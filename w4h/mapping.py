import datetime
import inspect
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

from w4h import logger_function

lidarURL = r'https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_Statewide_Lidar_DEM_WGS/ImageServer/WCSServer?request=GetCapabilities&service=WCS'

#Read study area shapefile (or other file) into geopandas
def read_study_area(study_area_path, study_area_crs='EPSG:4269', log=False):
    """Read study area geospatial file into geopandas

    Parameters
    ----------
    study_area_path : str or pathlib.Path
        Filepath to any geospatial file readable by geopandas. 
        Polygon is best, but may work with other types if extent is correct.
    crs : str, tuple, dict, optional
        CRS designation readable by geopandas/pyproj
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    studyAreaIN : geopandas dataframe
        Geopandas dataframe with polygon geometry.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    studyAreaIN = gpd.read_file(study_area_path)
    studyAreaIN.to_crs(study_area_crs, inplace=True)

    return studyAreaIN

#Convert coords in columns to geometry in geopandas dataframe
def coords2geometry(df_no_geometry, xcol='LONGITUDE', ycol='LATITUDE', zcol='ELEV_FT', input_coords_crs='EPSG:4269', use_z=False, log=False):
    """Adds geometry to points with xy coordinates in the specified coordinate reference system.

    Parameters
    ----------
    df_no_geometry : pandas.Dataframe
        a Pandas dataframe containing points
    xcol : str, default='LONGITUDE'
        Name of column holding x coordinate data in df_no_geometry
    ycol : str, default='LATITUDE'
        Name of column holding y coordinate data in df_no_geometry
    zcol : str, default='ELEV_FT'
        Name of column holding z coordinate data in df_no_geometry
    input_coords_crs : str, default='EPSG:4269
        Name of crs used for geometry
    use_z : bool, default=False
        Whether to use z column in calculation
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        Geopandas dataframe with points and their geometry values

    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    df = df_no_geometry.copy()

    ptCRS=input_coords_crs

    x = df[xcol].to_numpy()
    y = df[ycol].to_numpy()
    z = df[zcol].to_numpy()

    #coords = pd.concat([y, x], axis=1)
    if use_z:
        df["geometry"] = gpd.points_from_xy(x, y, z=z, crs=ptCRS)
    else:
        df["geometry"] = gpd.points_from_xy(x, y, crs=ptCRS)
        
    gdf = gpd.GeoDataFrame(df, crs=ptCRS)
    return gdf

#Clip a geodataframe to a study area
def clip_gdf2study_area(study_area, gdf, log=False):
    """Clips dataframe to only include things within study area.

    Parameters
    ----------
    study_area : geopandas.GeoDataFrame
        Inputs study area polygon
    gdf : geopandas.GeoDataFrame
        Inputs point data
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    gdfClip : geopandas.GeoDataFrame
        Contains only points within the study area
    
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if study_area is None:
        return gdf
    else:
        studyArea_proj = study_area.to_crs(gdf.crs).copy()
        gdfClip = gpd.clip(gdf, studyArea_proj) #Easier to project just study area to ensure data fit
        gdfClip.reset_index(inplace=True, drop=True) #Reset index
    
    return gdfClip

#Function to sample raster points to points specified in geodataframe
def sample_raster_points(raster, points_df, xcol='LONGITUDE', ycol='LATITUDE', new_col='SAMPLED', verbose=True, log=False):  
    """Sample raster values to points from geopandas geodataframe.

    Parameters
    ----------
    raster : rioxarray data array
        Raster containing values to be sampled.
    points_df : geopandas.geodataframe
        Geopandas dataframe with geometry column containing point values to sample.
    xcol : str, default='LONGITUDE'
        Column containing name for x-column, by default 'LONGITUDE.'
        This is used to output (potentially) reprojected point coordinates so as not to overwrite the original.
    ycol : str, default='LATITUDE'
        Column containing name for y-column, by default 'LATITUDE.'
        This is used to output (potentially) reprojected point coordinates so as not to overwrite the original.    new_col : str, optional
    new_col : str, default='SAMPLED'
        Name for name of new column containing points sampled from the raster, by default 'SAMPLED'.
    verbose : bool, default=True
        Whether to send to print() information about progress of function, by default True.
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    points_df : geopandas.geodataframe
        Same as points_df, but with sampled values and potentially with reprojected coordinates.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if verbose:
        nowTime = datetime.datetime.now()
        expectMin = (points_df.shape[0]/3054409) * 14
        endTime = nowTime+datetime.timedelta(minutes=expectMin)
        print(new_col+ " sampling should be done by {:d}:{:02d}".format(endTime.hour, endTime.minute))

    #Project points to raster CRS
    rastercrsWKT=raster.spatial_ref.crs_wkt
    points_df = points_df.to_crs(rastercrsWKT)
    #if xcol=='LONGITUDE' and ycol=='LATITUDE':
    xCOLOUT = xcol+'_PROJ'
    yCOLOUT = ycol+'_PROJ'
    points_df[xCOLOUT] = points_df['geometry'].x
    points_df[yCOLOUT] = points_df['geometry'].y
    xData = np.array(points_df[xCOLOUT].values)
    yData = np.array(points_df[yCOLOUT].values)
    sampleArr=raster.sel(x=xData, y=yData, method='nearest').values
    sampleArr = np.diag(sampleArr)
    sampleDF = pd.DataFrame(sampleArr, columns=[new_col])
    points_df[new_col] = sampleDF[new_col]
    return points_df

#Merge xyz dataframe into a metadata dataframe
def xyz_metadata_merge(xyz, metadata, log=False):
    """Add elevation to header data file.

    Parameters
    ----------
    xyz : pandas.Dataframe
        Contains elevation for the points
    metadata : pandas dataframe
        Header data file
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    headerXYZData : pandas.Dataframe
        Header dataset merged to get elevation values

    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    headerXYZData = metadata.merge(xyz, how='left', on='API_NUMBER')
    headerXYZData.drop(['LATITUDE_x', 'LONGITUDE_x'], axis=1, inplace=True)
    headerXYZData.rename({'LATITUDE_y':'LATITUDE', 'LONGITUDE_y':'LONGITUDE'}, axis=1, inplace=True)
    return headerXYZData

#Read wcsservice into rioxarray
def read_wcs(study_area, wcs_url=lidarURL, res_x=30, res_y=30, log=False, **kwargs):
    """Reads a WebCoverageService from a url and returns a rioxarray dataset containing it.

    Parameters
    ----------
    study_area : geopandas.GeoDataFrame
        Dataframe containing study area polygon
    wcs_url : str, default=lidarURL
    Represents the url for the WCS
    res_x : int, default=30
        Sets resolution for x axis
    res_y : int, default=30
        Sets resolution for y axis
    log : bool, default = False
        Whether to log results to log file, by default False
    **kwargs

    Returns
    -------
    wcsData_rxr : xarray.DataArray
        A xarray dataarray holding the image from the WebCoverageService
    """
    #Drawn largely from: https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/01-WCS-basics.ipynb
    
    #30m DEM
    #wcs_url = r'https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_DEM_30M/ImageServer/WCSServer?request=GetCapabilities&service=WCS'
    #lidar url:
    #lidarURL = r'https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_Statewide_Lidar_DEM_WGS/ImageServer/WCSServer?request=GetCapabilities&service=WCS'

    #studyAreaPath = r"\\isgs-sinkhole.ad.uillinois.edu\geophysics\Balikian\ISWS_HydroGeo\WellDataAutoClassification\SampleData\ESL_StudyArea_5mi.shp"
    #study_area = gpd.read_file(studyAreaPath)
    if study_area is None:
        print('ERROR: study_area must be specified to use read_wcs (currently set to {})'.format(study_area))
        return

    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if 'wcs_url' in kwargs:
        wcs_url = kwargs['wcs_url']
    if 'res_x' in kwargs:
        res_x = kwargs['res_x']
    if 'res_y' in kwargs:
        res_y = kwargs['res_y']
    
    width_in = ''
    height_in= ''

    #Create coverage object
    my_wcs = WebCoverageService(wcs_url, version='1.0.0') 
    #names = [k for k in my_wcs.contents.keys()]
    #print(names)
    dataID = 'IL_Statewide_Lidar_DEM'
    data = my_wcs.contents[dataID]
    dBBox = data.boundingboxes #Is this an error?
    
    study_area = study_area.to_crs(data.boundingboxes[0]['nativeSrs'])
    saBBox = study_area.total_bounds
    
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

#Read wms service into rioxarray
def read_wms(study_area, layer_name='IL_Statewide_Lidar_DEM_WGS:None', wms_url=lidarURL, srs='EPSG:3857', clip_to_studyarea=True, bbox=[-9889002.615500,5134541.069716,-9737541.607038,5239029.627400],res_x=30, res_y=30, size_x=512, size_y=512, format='image/tiff', log=False, **kwargs):
    """
    Reads a WebMapService from a url and returns a rioxarray dataset containing it.

    Parameters
    ----------
    study_area : geopandas.GeoDataFrame
        Dataframe containg study area polygon
    layer_name : str, default='IL_Statewide_Lidar_DEM_WGS:None'
        Represents the layer name in the WMS
    wms_url : str, default=lidarURL
        Represents the url for the WMS
    srs : str, default='EPSG:3857'
        Sets the srs
    clip_to_studyarea : bool, default=True
        Whether to clip to study area or not
    res_x : int, default=30
        Sets resolution for x axis
    res_y : int, default=512
        Sets resolution for y axis
    size_x : int, default=512
        Sets width of result
    size_y : int, default=512
        Sets height of result
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    wmsData_rxr : xarray.DataArray
        Holds the image from the WebMapService
    """
    if study_area is None:
        print('ERROR: study_area must be specified to use read_wms (currently set to {})'.format(study_area))
        return
    
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    from owslib.wms import WebMapService
    # Define WMS endpoint URL
    if 'wms_url' in kwargs:
        wms_url = kwargs['wms_url']
    else:
        wms_url = wms_url

    # Create WMS connection object
    wms = WebMapService(wms_url)
    # Print available layers
    #print(wms.contents)

    # Select desired layer
    if 'layer_name' in kwargs:
        layer = kwargs['layer_name']
    else:
        layer = layer_name
    
    data = wms.contents#[layer]
    if 'srs' in kwargs:
        studyArea_proj = study_area.to_crs(kwargs['srs'])
        saBBox = studyArea_proj.total_bounds
    else:
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
    # Check kwargs for rest of parameters
    if 'size_x' in kwargs:
        size_x = kwargs['size_x']
    if 'size_y' in kwargs:
        size_y = kwargs['size_y']
    if 'format' in kwargs:
        format = kwargs['format']
    if 'clip_to_studyarea' in kwargs:
        clip_to_studyarea = kwargs['clip_to_studyarea']
   
    #get the wms
    if clip_to_studyarea:
        img = wms.getmap(layers=[layer], srs=srs, bbox=saBBox, size=(size_x, size_y), format=format, transparent=True, timeout=60)        
    else:
        img = wms.getmap(layers=[layer], srs=srs, bbox=bbox, size=(size_x, size_y), format=format, transparent=True, timeout=60)

    #Save wms in memory to a raster dataset
    with MemoryFile(img) as memfile:
        with memfile.open() as dataset:
            wmsData_rxr = rxr.open_rasterio(dataset)

    #if clip_to_studyarea:
    #    wmsData_rxr = wmsData_rxr.sel(x=slice(saBBox[0], saBBox[2]), y=slice(saBBox[3], saBBox[1]))#.sel(band=1)

    return wmsData_rxr

#Clip a grid to a study area
def grid2study_area(study_area, grid, study_area_crs='', grid_crs='', log=False):
    """Clips grid to study area.

    Parameters
    ----------
    study_area : geopandas.GeoDataFrame
        inputs study area polygon
    grid : xarray.DataArray
        inputs grid array
    study_area_crs : str, default=''
        inputs the coordinate reference system for the study area
    grid_crs : str, default=''
        inputs the coordinate reference system for the grid
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    grid : xarray.DataArray
        returns xarray containing grid clipped only to area within study area

    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)
    
    if study_area_crs=='':
        study_area_crs=study_area.crs

    if grid_crs=='':
        #Get EPSG of model grid
        subtext = grid.spatial_ref.crs_wkt[-20:]
        starInd = subtext.find('EPSG')
        grid_crs = subtext[starInd:-2].replace('"','').replace(',',':')   
        #print(grid_crs)
    
    if study_area_crs != grid_crs:
        studyAreaUnproject = study_area.copy()
        study_area = study_area.to_crs(grid_crs)   
    else:
        study_area = study_area

    saExtent = study_area.total_bounds

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

#Read the model grid into (rio)xarray
def read_model_grid(model_grid_path, study_area=None, no_data_val_grid=0, read_grid=True, node_byspace=True, study_area_crs=None, grid_crs=None, verbose=False, log=False):
    """Reads in model grid to xarray data array

    Parameters
    ----------
    grid_path : str
        Path to model grid file
    study_area : geopandas.GeoDataFrame, default=None
        Dataframe containing study area polygon
    no_data_val_grid : int, default=0
        value assigned to areas with no data
    readGrid : bool, default=True
        Whether function to either read grid or create grid
    node_byspace : bool, default=False
        Denotes how to create grid
    study_area_crs : str, default=None
        Inputs study area crs
    grid_crs : str, default=None
        Inputs grid crs
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    modelGrid : xarray.DataArray
        Data array containing model grid
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)
    
    if read_grid and model_grid_path is not None:
        modelGridIN = rxr.open_rasterio(model_grid_path)

        file = w4h.read.__get_resource_path('isws_crs.txt')
        iswsCRS = w4h.read_dict(file, keytype=None)

        if grid_crs is None:
            try:
                grid_crs=modelGridIN.spatial_ref.crs_wkt
            except:
                modelGridIN.rio.write_crs(iswsCRS)
        elif grid_crs.lower()=='isws':
            modelGridIN.rio.write_crs(iswsCRS)
        
        if study_area is not None:                
            if study_area_crs is None:
                study_area_crs=study_area.crs
            study_area = study_area.to_crs(grid_crs)
            study_area_crs=study_area.crs            
            modelGrid = grid2study_area(study_area=study_area, grid=modelGridIN, study_area_crs=study_area_crs, grid_crs=grid_crs)
        else:
            modelGrid = modelGridIN

        try:
            noDataVal = float(modelGrid.attrs['_FillValue']) #Extract from dataset itsel
        except:
            noDataVal = -5000000

        modelGrid = modelGrid.where(modelGrid != noDataVal, other=np.nan)   #Replace no data values with NaNs
        modelGrid.rio.reproject(iswsCRS, inplace=True)
    elif model_grid_path is None and study_area is None:
        if verbose:
            print("ERROR: Either model_grid_path or study_area must be defined.")
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
        
        saExtent = study_area.total_bounds

        startX = saExtent[0] #Starting X Coordinate
        startY = saExtent[1] #starting Y Coordinate
        
        endX = saExtent[2]
        endY = saExtent[3]
        
        if node_byspace:
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
        if grid_crs is None or grid_crs=='isws' or grid_crs=='ISWS':
            for k in spatRefDict:
                modelGrid.spatial_ref.attrs[k] = spatRefDict[k]
    return modelGrid

#Read a grid from a file in using rioxarray
def read_grid(grid_path=None, grid_type='model', no_data_val_grid=0, use_service=False, study_area=None, study_area_crs=None, grid_crs=None, verbose=False, log=False, **kwargs):
    """Reads in grid

    Parameters
    ----------
    grid_path : str or pathlib.Path, default=None
        Path to a grid file
    grid_type : str, default='model'
        Sets what type of grid to load in
    no_data_val_grid : int, default=0
        Sets the no data value of the grid
    use_service : str, default=False
        Sets which service the function uses
    study_area : geopandas.GeoDataFrame, default=None
        Dataframe containing study area polygon
    study_area_crs : str, default=None
        Sets specific crs if current crs is not wanted
    grid_crs : str, default=None
        Sets crs to use if clipping to study area
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    gridIN : xarray.DataArray
        Returns grid
    
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if grid_type=='model':
        if 'read_grid' in list(kwargs.keys()):
            rgrid = kwargs['read_grid']
        else:
            rgrid=True
        gridIN = read_model_grid(model_grid_path=grid_path, study_area=study_area,  no_data_val_grid=0, read_grid=rgrid, study_area_crs=study_area_crs, grid_crs=grid_crs, verbose=verbose)
    else:
        if use_service==False:
            gridIN = rxr.open_rasterio(grid_path)
        elif use_service.lower()=='wcs':
            gridIN = read_wcs(study_area, wcs_url=lidarURL, **kwargs)
        elif use_service.lower()=='wms':
            pass
            gridIN = read_wms(study_area, wcs_url=lidarURL, **kwargs)
            
        if study_area is not None:
            if grid_crs is None:
                try:
                    grid_crs=gridIN.spatial_ref.crs_wkt
                except:
                    iswsCRS = w4h.read_dict(r'../resources/isws_crs')
                    gridIN.rio.write_crs(iswsCRS)
            elif grid_crs.lower()=='isws':
                iswsCRS = w4h.read_dict(r'../resources/isws_crs')
                gridIN.rio.write_crs(iswsCRS)
                
            if study_area_crs is None:
                study_area_crs=study_area.crs
            study_area = study_area.to_crs(grid_crs)
            study_area_crs=study_area.crs
            
            gridIN = grid2study_area(study_area=study_area, grid=gridIN, study_area_crs=study_area_crs, grid_crs=grid_crs)

        try:
            no_data_val_grid = gridIN.attrs['_FillValue'] #Extract from dataset itself
        except:
            pass
                
        gridIN = gridIN.where(gridIN != no_data_val_grid, other=np.nan)  #Replace no data values with NaNs

    return gridIN

#Align and coregister rasters
def align_rasters(grids_unaligned, modelgrid, no_data_val_grid=0, log=False):
    """Reprojects two rasters and aligns their pixels

    Parameters
    ----------
    grids_unaligned : list or xarray.DataArray
        Contains a list of grids or one unaligned grid
    modelgrid : xarray.DataArray
        Contains model grid
    no_data_val_grid : int, default=0
        Sets value of no data pixels
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    alignedGrids : list or xarray.DataArray
        Contains aligned grids
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if type(grids_unaligned) is list:
        alignedGrids=[]
        for g in grids_unaligned:
            alignedGrid = g.rio.reproject_match(modelgrid)

            try:
                no_data_val_grid = alignedGrid.attrs['_FillValue'] #Extract from dataset itself
            except:
                pass
            
            alignedGrid = alignedGrid.where(alignedGrid != no_data_val_grid)  #Replace no data values with NaNs
            
            alignedGrids.append(alignedGrid)
    else:
        alignedGrid = grids_unaligned.rio.reproject_match(modelgrid)

        try:
            noDataVal = alignedGrid.attrs['_FillValue'] #Extract from dataset itself
        except:
            pass

        alignedGrids = alignedGrid.where(alignedGrid != noDataVal, other=np.nan)  #Replace no data values with NaNs
        
    return alignedGrids

#Get drift and layer thickness, given a surface and bedrock grid
def get_drift_thick(surface, bedrock, layers=9, plot=False, log=False):
    """Finds the distance from surface to bedrock and then divides by number of layers to get layer thickness.

    Parameters
    ----------
    surface : rioxarray.DataArray
        array holding surface elevation
    bedrock : rioxarray.DataArray
        array holding bedrock elevation
    layers : int, default=9
        number of layers needed to calculate thickness for
    plot : bool, default=False
        tells function to either plot the data or not

    Returns
    -------
    driftThick : rioxarray.DataArray
        Contains data array containing depth to bedrock at each point
    layerThick : rioxarray.DataArray
        Contains data array with layer thickness at each point

    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    xr.set_options(keep_attrs=True)

    driftThick = surface - bedrock
    driftThick = driftThick.clip(0,max=5000,keep_attrs=True)
    if plot:
        driftThick.plot()

    try:
        noDataVal = driftThick.attrs['_FillValue'] #Extract from dataset itself
    except:
        noDataVal = 100001
    
    driftThick = driftThick.where(driftThick <100000, other=np.nan)  #Replace no data values with NaNs
    driftThick = driftThick.where(driftThick >-100000, other=np.nan)  #Replace no data values with NaNs

    layerThick = driftThick/layers
    
    xr.set_options(keep_attrs='default')

    return driftThick, layerThick