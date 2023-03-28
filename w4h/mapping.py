import datetime
import json
import pathlib
import os


import rioxarray as rxr
import xarray as xr
import geopandas as gpd
import pandas as pd
import numpy as np
import shapely
from owslib.wcs import WebCoverageService
from shapely.geometry import Point
from urllib.request import urlopen
from rasterio import MemoryFile

import w4h

lidarURL = r'https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_Statewide_Lidar_DEM_WGS/ImageServer/WCSServer?request=GetCapabilities&service=WCS'

#Read study area shapefile (or other file) into geopandas
def read_study_area(studyareapath, crs=''):
    """Read study area geospatial file into geopandas

    Parameters
    ----------
    studyareapath : str or pathlib.Path
        Filepath to any geospatial file readable by geopandas. 
        Polygon is best, but may work with other types if extent is correct.
    crs : str, tuple, dict, optional
        CRS designation readable by geopandas/pyproj

    Returns
    -------
    studyAreaIN : geopandas dataframe
        Geopandas dataframe with polygon geometry.
    """
    studyAreaIN = gpd.read_file(studyareapath)
    return studyAreaIN


def coords2Geometry(df, xCol='LONGITUDE', yCol='LATITUDE', zCol='ELEV_FT', crs='EPSG:4269', useZ=False):
    '''
    
    Adds geometry to points with xy coordinates in the specified coordinate reference system.

            Parameters:
                    df (pandas dataframe): a Pandas dataframe containing points
                    xCol (str): Name of column holding x coordinate data in df
                    yCol (str): Name of column holding y coordinate data in df
                    zCol (str): Name of column holding z coordinate data in df
                    crs (str): Name of crs used for geometry
                    useZ (bool): Whether to use z column in calculation

            Returns:
                    gdf (geopandas dataframe): Geopandas dataframe with points and their geometry values

    '''

    ptCRS=crs

    x = df[xCol].to_numpy()
    y = df[yCol].to_numpy()
    z = df[zCol].to_numpy()

    #coords = pd.concat([y, x], axis=1)
    if useZ:
        df["geometry"] = gpd.points_from_xy(x, y, z=z, crs=ptCRS)
    else:
        df["geometry"] = gpd.points_from_xy(x, y, crs=ptCRS)
        
    gdf = gpd.GeoDataFrame(df, crs=ptCRS)
    return gdf

def clipHeader2StudyArea(studyarea, headerdata, headerCRS='EPSG:4269'):
    '''
    
    Clips dataframe to only include things within study area.

            Parameters:
                    studyarea (geopandas dataframe): Inputs study area polygon
                    headerdata (geopandas dataframe): Inputs point data
                    headerCRS (str): Inputs crs to project study area to
            
            Returns:
                    headerDataClip (geopandas dataframe): Contains only points within the study area
    
    '''
    studyArea_4269 = studyarea.to_crs(headerCRS).copy()
    
    headerDataClip = gpd.clip(headerdata, studyArea_4269) #Data from table is in EPSG:4269, easier to just project study area to ensure data fit
    
    headerDataClip.reset_index(inplace=True, drop=True) #Reset index
    
    return headerDataClip

def sample_raster_points(raster, ptDF, xCol='LONGITUDE', yCol='LATITUDE', newColName='SAMPLED', printouts=True):  
    """Sample raster values to points from geopandas geodataframe.

    Parameters
    ----------
    raster : rioxarray data array
        Raster containing values to be sampled.
    ptDF : geopandas.geodataframe
        Geopandas dataframe with geometry column containing point values to sample.
    xCol : str, default='LONGITUDE'
        Column containing name for x-column, by default 'LONGITUDE.'
        This is used to output (potentially) reprojected point coordinates so as not to overwrite the original.
    yCol : str, default='LATITUDE'
        Column containing name for y-column, by default 'LATITUDE.'
        This is used to output (potentially) reprojected point coordinates so as not to overwrite the original.    newColName : str, optional
    newColName : str, default='SAMPLED'
        Name for name of new column containing points sampled from the raster, by default 'SAMPLED'.
    printouts : bool, default=True
        Whether to send to print() information about progress of function, by default True.

    Returns
    -------
    ptDF : geopandas.geodataframe
        Same as ptDF, but with sampled values and potentially with reprojected coordinates.
    """
    if printouts:
        nowTime = datetime.datetime.now()
        expectMin = (ptDF.shape[0]/3054409) * 14
        endTime = nowTime+datetime.timedelta(minutes=expectMin)
        print(newColName+ " sampling should be done by {:d}:{:02d}".format(endTime.hour, endTime.minute))

    #Project points to raster CRS
    rastercrsWKT=raster.spatial_ref.crs_wkt
    ptDF = ptDF.to_crs(rastercrsWKT)
    #if xCol=='LONGITUDE' and yCol=='LATITUDE':
    xCOLOUT = xCol+'_PROJ'
    yCOLOUT = yCol+'_PROJ'
    ptDF[xCOLOUT] = ptDF['geometry'].x
    ptDF[yCOLOUT] = ptDF['geometry'].y
    xData = np.array(ptDF[xCOLOUT].values)
    yData = np.array(ptDF[yCOLOUT].values)
    sampleArr=raster.sel(x=xData, y=yData, method='nearest').values
    sampleArr = np.diag(sampleArr)
    sampleDF = pd.DataFrame(sampleArr, columns=[newColName])
    ptDF[newColName] = sampleDF[newColName]
    return ptDF

def addElevtoHeader(xyz, header):
    '''
    
    Adds elevation to header data file.

            Parameters:
                    xyz (pandas dataframe): Contains elevation for the points
                    header (pandas dataframe): Header data file
                
            Returns:
                    headerXYZData (pandas dataframe): Header dataset merged to get elevation values
    
    '''
    headerXYZData = header.merge(xyz, how='left', on='API_NUMBER')
    headerXYZData.drop(['LATITUDE_x', 'LONGITUDE_x'], axis=1, inplace=True)
    headerXYZData.rename({'LATITUDE_y':'LATITUDE', 'LONGITUDE_y':'LONGITUDE'}, axis=1, inplace=True)
    return headerXYZData

def readWCS(studyArea, wcs_url=lidarURL, res_x=30, res_y=30):
    '''
    
    Reads a WebCoverageService from a url and returns a rioxarray dataset containing it.

            Parameters:
                    studyArea (list): A list of integers representing a bounding box
                    wcs_url (str): A string representing the url for the WCS
                    res_x (int): An integer to set resolution for x axis
                    res_y (int): An integer to set resolution for y axis
            
            Returns:
                    wcsData_rxr (xarray dataarray): A xarray dataarray holding the image from the WebCoverageService


    '''
    #Drawn largely from: https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/01-WCS-basics.ipynb
    
    #30m DEM
    #wcs_url = r'https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_DEM_30M/ImageServer/WCSServer?request=GetCapabilities&service=WCS'
    #lidar url:
    #lidarURL = r'https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_Statewide_Lidar_DEM_WGS/ImageServer/WCSServer?request=GetCapabilities&service=WCS'

    #studyAreaPath = r"\\isgs-sinkhole.ad.uillinois.edu\geophysics\Balikian\ISWS_HydroGeo\WellDataAutoClassification\SampleData\ESL_StudyArea_5mi.shp"
    #studyArea = gpd.read_file(studyAreaPath)
    
    width_in = ''
    height_in= ''

    #Create coverage object
    my_wcs = WebCoverageService(wcs_url, version='1.0.0') 
    #names = [k for k in my_wcs.contents.keys()]
    #print(names)
    dataID = 'IL_Statewide_Lidar_DEM'
    data = my_wcs.contents[dataID]
    dBBox = data.boundingboxes #Is this an error?
    
    studyArea = studyArea.to_crs(data.boundingboxes[0]['nativeSrs'])
    saBBox = studyArea.total_bounds
    
    #In case study area bounding box goes outside data bounding box, use data bounding box values
    newBBox = []
    for i,c in enumerate(dBBox[0]['bbox']):
        if i == 0 or i==2:
            if saBBox[i] < c:
                newBBox.append(saBBox[i])
            else:
                newBBox.append(c)
        else:
            if saBBox[i] > c:
                newBBox.append(saBBox[i])
            else:
                newBBox.append(c)

    #Recalculate resolution if it is too fine to read in
    #Start by getting the area of the study area bounding box
    saWidth = saBBox[2]-saBBox[0]
    saHeight = saBBox[3]-saBBox[1]
    saBBoxAreaM = saWidth*saHeight
    saBBoxAreaKM = saBBoxAreaM/(1000*1000) #Area in km^2

    if saBBoxAreaM/(res_x*res_y) > (4100*15000)*0.457194: #What I think might be the max download size?
        print("Resolution inputs overriden, file request too large.")
        res_x=str(round(saWidth/2500, 2))

        width_in  = str(int(saWidth/float(res_x )))
        height_in = str(int(saHeight/float(res_x)))
        
        res_y=str(round(saHeight/height_in, 2))

        print('New resolution is: '+res_x+'m_x X '+res_y+'m_y' )
        print('Dataset size: '+width_in+' pixels_x X '+height_in+' pixels_y')

    bBox = tuple(newBBox)
    bBox_str = str(tuple(newBBox)[1:-1]).replace(' ','')
    dataCRS = 'EPSG:3857'

    #Format WCS request using owslib
    response = my_wcs.getCoverage(
        identifier=my_wcs[dataID].id, 
        crs=dataCRS,#'urn:ogc:def:crs:EPSG::26716',
        bbox=bBox,
        resx=res_x, 
        resy=res_y,
        timeout=60,
        #width = width_in, height=height_in,
        format='GeoTIFF')
    response

    #If I can figure out url, this might be better?
    #baseURL = r'https://data.isgs.illinois.edu/arcgis/services/Elevation/'+dataID+'/ImageServer/WCSServer'
    #addonRequestURL = '?request=GetCoverage&service=WCS&bbox='+bBox_str+'&srs='+dataCRS+'&format=GeoTIFF'+'&WIDTH='+width_in+'&HEIGHT='+height_in+')'
    #reqURL = baseURL+addonRequestURL
    #wcsData_rxr =  rxr.open_rasterio(reqURL)

    with MemoryFile(response) as memfile:
        with memfile.open() as dataset:
            wcsData_rxr =  rxr.open_rasterio(dataset)

    return wcsData_rxr

def readWMS(study_area, layer_name='IL_Statewide_Lidar_DEM_WGS:None', wms_url=lidarURL, srs='EPSG:3857', clip_to_studyarea=True, bbox=[-9889002.615500,5134541.069716,-9737541.607038,5239029.627400],size_x=256, size_y=256, format='image/tiff'):
    '''
    Reads a WebMapService from a url and returns a rioxarray dataset containing it.

            Parameters:
                    study_area (list): A list of integers representing a bounding box
                    layer_name (str): A string representing the layer name in the WMS
                    wms_url (str): A string representing the url for the WMS
                    res_x (int): An integer to set resolution for x axis
                    res_y (int): An integer to set resolution for y axis
                    size_x (int): An integer to set width of result
                    size_y (int): An integer to set height of result
            
            Returns:
                    wmsData_rxr (xarray dataarray): A xarray dataarray holding the image from the WebMapService
    '''
    from owslib.wms import WebMapService
    # Define WMS endpoint URL
    wms_url = wms_url
    # Create WMS connection object
    wms = WebMapService(wms_url)
    # Print available layers
    #print(wms.contents)
    # Select desired layer 
    layer = layer_name
    
    data = wms.contents#[layer]
    #print(data['0'].__dict__)
    studyArea_proj = study_area.to_crs(srs)
    saBBox = studyArea_proj.total_bounds
    
    if layer == 'IL_Statewide_Lidar_DEM_WGS:None':
        dBBox = data['0'].boundingBox #Is this an error?

        gpdDict = {'Label': ['Surf Data Box'], 'geometry': [shapely.geometry.Polygon(((dBBox[0], dBBox[1]), (dBBox[0], dBBox[3]), (dBBox[2], dBBox[3]), (dBBox[2], dBBox[1]), (dBBox[0], dBBox[1])))]}
        dBBoxGDF = gpd.GeoDataFrame(gpdDict, crs=dBBox[4])
        dBBoxGDF.to_crs(srs)
        #In case study area bounding box goes outside data bounding box, use data bounding box values
        newBBox = []
        for i,c in enumerate(dBBox):
            if type(c) is str:
                pass
            elif i == 0 or i==2:
                if saBBox[i] < c:
                    newBBox.append(saBBox[i])
                else:
                    newBBox.append(c)
            else:
                if saBBox[i] > c:
                    newBBox.append(saBBox[i])
                else:
                    newBBox.append(c)

    saWidth = saBBox[2]-saBBox[0]
    saHeight = saBBox[3]-saBBox[1]    
    #get the wms
    if clip_to_studyarea:
        img = wms.getmap(layers=[layer], srs=srs, bbox=saBBox, size=(size_x, size_y), format=format, transparent=True)        
    else:
        img = wms.getmap(layers=[layer], srs=srs, bbox=bbox, size=(size_x, size_y), format=format, transparent=True)


    #Save wms in memory to a raster dataset
    with MemoryFile(img) as memfile:
        with memfile.open() as dataset:
            wmsData_rxr = rxr.open_rasterio(dataset)
    wmsData_rxr = wmsData_rxr.astype(np.half)

    return wmsData_rxr

def clipGrid2StudyArea(studyArea, grid, studyAreacrs='', gridcrs=''):
    '''
    
    Clips grid to study area.

            Parameters:
                    studyArea (geopandas dataframe): inputs study area polygon
                    grid (xarray dataarray): inputs grid array
                    studyAreacrs (str): inputs the coordinate reference system for the study area
                    gridcrs (str): inputs the coordinate reference system for the grid
            
            Returns:
                    grid (xarray dataarray): returns xarray containing grid clipped only to area within study area

    '''
    
    if studyAreacrs=='':
        studyAreacrs=studyArea.crs

    if gridcrs=='':
        #Get EPSG of model grid
        subtext = grid.spatial_ref.crs_wkt[-20:]
        starInd = subtext.find('EPSG')
        gridcrs = subtext[starInd:-2].replace('"','').replace(',',':')   
        #print(gridcrs)
    
    if studyAreacrs != gridcrs:
        studyAreaUnproject = studyArea.copy()
        studyArea = studyArea.to_crs(gridcrs)   
    else:
        studyArea = studyArea

    saExtent = studyArea.total_bounds

    if grid['y'][-1].values - grid['y'][0].values > 0:
        miny=saExtent[1]
        maxy=saExtent[3]
    else:
        miny=saExtent[3]
        maxy=saExtent[1]        
        
    if grid['x'][-1].values - grid['x'][0].values > 0:
        minx=saExtent[0]
        maxx=saExtent[2]
    else:
        minx=saExtent[2]
        maxx=saExtent[0]
    grid = grid.sel(x=slice(minx, maxx), y=slice(miny, maxy)).sel(band=1)     

    return grid


def read_model_grid(studyArea, gridpath, nodataval=0, read_grid=True, node_bySpace=False, clip2SA=True, studyAreacrs=None, gridcrs=None):
    """_summary_

    Parameters
    ----------
    studyArea : _type_
        _description_
    gridpath : _type_
        _description_
    nodataval : int, optional
        _description_, by default 0
    readGrid : bool, optional
        _description_, by default True
    node_bySpace : bool, optional
        _description_, by default False
    clip2SA : bool, optional
        _description_, by default True
    studyAreacrs : str, optional
        _description_, by default ''
    gridcrs : str, optional
        _description_, by default ''

    Returns
    -------
    _type_
        _description_
    """
    read_grid = True
    node_bySpace = False #False means by number of nodes
    
    if read_grid:
        modelGridIN = rxr.open_rasterio(gridpath)

        file = w4h.get_resource_path('isws_crs.txt')
        iswsCRS = w4h.read_dict(file, keytype=None)

        if gridcrs is None:
            try:
                gridcrs=modelGridIN.spatial_ref.crs_wkt
            except:
                modelGridIN.rio.write_crs(iswsCRS)
        elif gridcrs.lower()=='isws':
            modelGridIN.rio.write_crs(iswsCRS)
        
        if clip2SA:                
            if studyAreacrs is None:
                studyAreacrs=studyArea.crs
            studyArea = studyArea.to_crs(gridcrs)
            studyAreacrs=studyArea.crs            
            modelGrid = clipGrid2StudyArea(studyArea=studyArea, grid=modelGridIN, studyAreacrs=studyAreacrs, gridcrs=gridcrs)
        try:
            noDataVal = float(modelGrid.attrs['_FillValue']) #Extract from dataset itsel
        except:
            noDataVal = -5000000

        modelGrid = modelGrid.where(modelGrid != noDataVal, other=np.nan)   #Replace no data values with NaNs
        modelGrid.rio.reproject(iswsCRS, inplace=True)
    else:
        spatRefDict = {'crs_wkt': 'PROJCS["Clarke_1866_Lambert_Conformal_Conic",GEOGCS["NAD27",DATUM["North_American_Datum_1927",SPHEROID["Clarke 1866",6378206.4,294.978698199999,AUTHORITY["EPSG","7008"]],AUTHORITY["EPSG","6267"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4267"]],PROJECTION["Lambert_Conformal_Conic_2SP"],PARAMETER["latitude_of_origin",33],PARAMETER["central_meridian",-89.5],PARAMETER["standard_parallel_1",33],PARAMETER["standard_parallel_2",45],PARAMETER["false_easting",2999994],PARAMETER["false_northing",0],UNIT["US survey foot",0.304800609601219,AUTHORITY["EPSG","9003"]],AXIS["Easting",EAST],AXIS["Northing",NORTH]]',
            'semi_major_axis': 6378206.4,
            'semi_minor_axis': 6356583.799998981,
            'inverse_flattening': 294.978698199999,
            'reference_ellipsoid_name': 'Clarke 1866',
            'longitude_of_prime_meridian': 0.0,
            'prime_meridian_name': 'Greenwich',
            'geographic_crs_name': 'NAD27',
            'horizontal_datum_name': 'North American Datum 1927',
            'projected_crs_name': 'Clarke_1866_Lambert_Conformal_Conic',
            'grid_mapping_name': 'lambert_conformal_conic',
            'standard_parallel': (33.0, 45.0),
            'latitude_of_projection_origin': 33.0,
            'longitude_of_central_meridian': -89.5,
            'false_easting': 2999994.0,
            'false_northing': 0.0,
            'spatial_ref': 'PROJCS["Clarke_1866_Lambert_Conformal_Conic",GEOGCS["NAD27",DATUM["North_American_Datum_1927",SPHEROID["Clarke 1866",6378206.4,294.978698199999,AUTHORITY["EPSG","7008"]],AUTHORITY["EPSG","6267"]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4267"]],PROJECTION["Lambert_Conformal_Conic_2SP"],PARAMETER["latitude_of_origin",33],PARAMETER["central_meridian",-89.5],PARAMETER["standard_parallel_1",33],PARAMETER["standard_parallel_2",45],PARAMETER["false_easting",2999994],PARAMETER["false_northing",0],UNIT["US survey foot",0.304800609601219,AUTHORITY["EPSG","9003"]],AXIS["Easting",EAST],AXIS["Northing",NORTH]]',
            'GeoTransform': '2440250.0 625.0 0.0 3459750.0 0.0 -625.0'}
        
        saExtent = studyArea.total_bounds

        startX = saExtent[0] #Starting X Coordinate
        startY = saExtent[1] #starting Y Coordinate
        
        endX = saExtent[2]
        endY = saExtent[3]
        
        if node_bySpace:
            xSpacing = 625 #X Node spacing 
            ySpacing = xSpacing #Y Node spacing  
            
            x = np.arange(startX, endX, xSpacing)
            y = np.arange(startY, endY, ySpacing)
        else:
            xNodes = 100 #Number of X Nodes
            yNodes = 100 #Number of Y Nodes

            x = np.linspace(startX, endX, num=xNodes)
            y = np.linspace(startY, endY, num=yNodes)        
        
        xx, yy = np.meshgrid(x, y)
        zz = np.ones_like(xx).transpose()

        yIn = np.flipud(y)

        coords = {'x':x,'y':yIn, 'spatial_ref':0}
        dims = {'x':x,'y':yIn}
        
        modelGrid = xr.DataArray(data=zz,coords=coords,attrs={'_FillValue':3.402823466e+38}, dims=dims)
        modelGrid.spatial_ref.attrs['spatial_ref'] = {}
        if gridcrs is None or gridcrs=='isws' or gridcrs=='ISWS':
            for k in spatRefDict:
                modelGrid.spatial_ref.attrs[k] = spatRefDict[k]
    return modelGrid


def read_grid(datapath, grid_type='model', nodataval=0, use_service=None, studyArea='', clip2SA=True,  studyAreacrs=None, gridcrs=None, **kwargs):
    if grid_type=='model':
        if 'read_grid' in list(kwargs.keys()):
            rgrid = kwargs['read_grid']
        else:
            rgrid=True
        gridIN = read_model_grid(studyArea, gridpath=datapath, nodataval=0, read_grid=rgrid, clip2SA=clip2SA, studyAreacrs=studyAreacrs, gridcrs=gridcrs)
    else:
        if use_service is None:
            gridIN = rxr.open_rasterio(datapath)
        elif use_service.lower()=='wcs':
            gridIN = readWCS(studyArea, wcs_url=lidarURL, **kwargs)
        elif use_service.lower()=='wms':
            gridIN = readWMS(studyArea, wms_url=lidarURL, **kwargs)
            
        if clip2SA:
            if gridcrs is None:
                try:
                    gridcrs=gridIN.spatial_ref.crs_wkt
                except:
                    iswsCRS = w4h.read_dict(r'../resources/isws_crs')
                    gridIN.rio.write_crs(iswsCRS)
            elif gridcrs.lower()=='isws':
                iswsCRS = w4h.read_dict(r'../resources/isws_crs')
                gridIN.rio.write_crs(iswsCRS)
                
            if studyAreacrs is None:
                studyAreacrs=studyArea.crs
            studyArea = studyArea.to_crs(gridcrs)
            studyAreacrs=studyArea.crs
            
            gridIN = clipGrid2StudyArea(studyArea=studyArea, grid=gridIN, studyAreacrs=studyAreacrs, gridcrs=gridcrs)
        try:
            nodataval = gridIN.attrs['_FillValue'] #Extract from dataset itself
        except:
            pass
                
        gridIN = gridIN.where(gridIN != nodataval, other=np.nan)  #Replace no data values with NaNs

    return gridIN

def alignRasters(unalignedGrids, modelgrid, nodataval=0):
    
    if type(unalignedGrids) is list:
        alignedGrids=[]
        for g in unalignedGrids:
            alignedGrid = g.rio.reproject_match(modelgrid)

            try:
                nodataval = alignedGrid.attrs['_FillValue'] #Extract from dataset itself
            except:
                pass
            
            alignedGrid = alignedGrid.where(alignedGrid != nodataval)  #Replace no data values with NaNs
            
            alignedGrids.append(alignedGrid)
    else:
        alignedGrid = unalignedGrids.rio.reproject_match(modelgrid)

        try:
            noDataVal = alignedGrid.attrs['_FillValue'] #Extract from dataset itself
        except:
            pass

        alignedGrids = alignedGrid.where(alignedGrid != noDataVal, other=np.nan)  #Replace no data values with NaNs
        
    return alignedGrids

def get_drift_thick(surface, bedrock, noLayers=9, plotData=False):
    '''
    
    Finds the distance from surface to bedrock and then divides by number of layers to get layer thickness.

            Parameters:
                    surface (rioxarray dataarray): array holding surface elevation
                    bedrock (rioxarray dataarray): array holding bedrock elevation
                    noLayers (int): number of layers needed to calculate thickness for
                    plotData (bool): tells function to either plot the data or not

            Returns:
                    driftThick (rioxarray dataarray): Contains data array containing depth to bedrock at each point
                    layerThick (rioxarray dataarray): Contains data array with layer thickness at each point
    
    '''
    xr.set_options(keep_attrs=True)

    driftThick = surface - bedrock
    driftThick = driftThick.clip(0,max=5000,keep_attrs=True)
    if plotData:
        driftThick.plot()

    try:
        noDataVal = driftThick.attrs['_FillValue'] #Extract from dataset itself
    except:
        noDataVal = 100001
    


    driftThick = driftThick.where(driftThick <100000, other=np.nan)  #Replace no data values with NaNs
    driftThick = driftThick.where(driftThick >-100000, other=np.nan)  #Replace no data values with NaNs

    layerThick = driftThick/noLayers
    
    xr.set_options(keep_attrs='default')

    return driftThick, layerThick
    