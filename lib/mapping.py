import rioxarray as rxr
import xarray as xr
import geopandas as gpd
import pandas as pd
import numpy as np
from owslib.wcs import WebCoverageService
from shapely.geometry import Point
import datetime
from urllib.request import urlopen
from rasterio import MemoryFile

lidarURL = r'https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_Statewide_Lidar_DEM_WGS/ImageServer/WCSServer?request=GetCapabilities&service=WCS'


def readStudyArea(studyareapath):
    studyAreaIN = gpd.read_file(studyareapath)
    saExtent = studyAreaIN.total_bounds
    return studyAreaIN, saExtent

def coords2Geometry(df, xCol='LONGITUDE', yCol='LATITUDE', zCol='ELEV_FT', crs='EPSG:4269', useZ=False):
    ptCRS=crs

    y = df[yCol]
    x = df[xCol]
    z = df[zCol]

    #coords = pd.concat([y, x], axis=1)
    if useZ:
        df["geometry"] = gpd.points_from_xy(y, x, z=z, crs=ptCRS)
    else:
        df["geometry"] = gpd.points_from_xy(y, x, crs=ptCRS)
        
    gdf = gpd.GeoDataFrame(df, crs=ptCRS)
    return gdf

def rastertoPoints_sample(raster, ptDF, xCol='LONGITUDE', yCol='LATITUDE', newColName='Sampled', printouts=True):
    if printouts:
        nowTime = datetime.datetime.now()
        endTime = nowTime+datetime.timedelta(minutes=14)
        print("Sampling process should be done by {:d}:{:02d}".format(endTime.hour, endTime.minute))

    if 'geometry' in ptDF.columns:
        if xCol=='LONGITUDE' and yCol=='LATITUDE':
            ptDF['LONGITUDE_PROJ'] = ptDF['geometry'].x
            ptDF['LATITUDE_PROJ'] = ptDF['geometry'].y
            xData = ptDF['LONGITUDE_PROJ']
            yData = ptDF['LATITUDE_PROJ']
    else:
        xData = ptDF['LONGITUDE']
        yData = ptDF['LATITUDE']

    sampleList = []
    for p in range(ptDF.shape[0]):
        sampleList.append(raster.sel(x=xData[p], y=yData[p], method='nearest').values[0])

    sampleDF = pd.DataFrame(sampleList, columns=[newColName])
    ptDF[newColName] = sampleDF[newColName]
    return ptDF

def addElevtoHeader(xyz, header):
    headerXYZData = header.merge(xyz, how='left', on='API_NUMBER')
    headerXYZData.drop(['LATITUDE_x', 'LONGITUDE_x'], axis=1, inplace=True)
    headerXYZData.rename({'LATITUDE_y':'LATITUDE', 'LONGITUDE_y':'LONGITUDE'}, axis=1, inplace=True)
    return headerXYZData

def readWCS(studyArea, wcs_url=lidarURL, res_x=30, res_y=30):

    #30m DEM
    #wcs_url = r'https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_DEM_30M/ImageServer/WCSServer?request=GetCapabilities&service=WCS'
    #lidar url:
    #lidarURL = r'https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_Statewide_Lidar_DEM_WGS/ImageServer/WCSServer?request=GetCapabilities&service=WCS'

    #studyAreaPath = r"\\isgs-sinkhole.ad.uillinois.edu\geophysics\Balikian\ISWS_HydroGeo\WellDataAutoClassification\SampleData\ESL_StudyArea_5mi.shp"
    #studyArea = gpd.read_file(studyAreaPath)
    
    studyArea = studyArea.to_crs(data.boundingboxes[0]['nativeSrs'])
    saBBox = studyArea.total_bounds

    width_in = ''
    height_in= ''

    #Create coverage object
    my_wcs = WebCoverageService(wcs_url, version='1.0.0') 
    #names = [k for k in my_wcs.contents.keys()]
    #print(names)
    dataID = 'IL_Statewide_Lidar_DEM'
    data = my_wcs.contents[dataID]
    dBBox = data.boundingboxes

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
    #wcsData_rxr =  rxr.open_rasterio(dataset)


    with MemoryFile(response) as memfile:
        with memfile.open() as dataset:
            wcsData_rxr =  rxr.open_rasterio(dataset)

    return wcsData_rxr

def readModelGrid(studyArea, gridpath, nodataval=0, readGrid=True, node_bySpace=False, clip2SA=True):
    
    saExtent = studyArea.total_bounds
    readGrid = True
    node_bySpace = False #False means by number of nodes

    if readGrid:
        modelGridIN = rxr.open_rasterio(gridpath)
        if clip2SA:
            modelGrid = modelGridIN.sel(x=slice(saExtent[0], saExtent[2]), y=slice(saExtent[1], saExtent[3])).sel(band=1)
        noDataMod = modelGrid.attrs['_FillValue'] #Extract from dataset itself
        modelGrid = modelGrid.where(modelGrid > noDataMod)   #Replace no data values with NaNs
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
        for k in spatRefDict:
            modelGrid.spatial_ref.attrs[k] = spatRefDict[k]

    return modelGrid

def readSurfaceGrid(surfaceelevpath='', nodataval=0, useWCS=False, studyArea=''):
    if useWCS:
        surfaceElevGridIN = readWCS(studyArea, wcs_url=lidarURL)
    else:
        surfaceElevGridIN = rxr.open_rasterio(surfaceelevpath)
        
        try:
            noDataSurf = surfaceElevGridIN.attrs['_FillValue'] #Extract from dataset itself
        except:
            if noDataSurf < -1000:
                noDataSurf = -1000
            elif noDataSurf > 50000:
                noDataSurf = 50000
            else:
                noDataSurf = nodataval #apply no data value
            
        surfaceElevGridIN = surfaceElevGridIN.where(surfaceElevGridIN != noDataSurf)  #Replace no data values with NaNs
        
    return surfaceElevGridIN

def readBedrockGrid(bedrockelevpath, nodataval=0):
    bedrockElevGridIN = rxr.open_rasterio(bedrockelevpath)
    
    try:
        noDataBed = bedrockElevGridIN.attrs['_FillValue'] #Extract from dataset itself
    except:
        if noDataBed < -1000:
            noDataBed = -1000
        elif noDataBed > 50000:
            noDataBed = 50000
        else:
            noDataBed = nodataval #apply no data value

    bedrockElevGridIN = bedrockElevGridIN.where(bedrockElevGridIN != noDataBed)   #Replace no data values with NaNs
    return bedrockElevGridIN

def alignRasters(bedrockgrid, surfacegrid, modelgrid, nodataval=0):
    surfaceGridAlign = surfacegrid.rio.reproject_match(modelgrid)
    bedrockGridAlign = bedrockgrid.rio.reproject_match(modelgrid)

    try:
        noDataBed = bedrockGridAlign.attrs['_FillValue'] #Extract from dataset itself
    except:
        if noDataBed < -1000:
            noDataBed = -1000
        elif noDataBed > 50000:
            noDataBed = 50000
        else:
            noDataBed = nodataval #apply no data value    
            
    try:
        noDataBed = surfaceGridAlign.attrs['_FillValue'] #Extract from dataset itself
    except:
        if noDataSurf < -1000:
            noDataSurf = -1000
        elif noDataSurf > 50000:
            noDataSurf = 50000
        else:
            noDataSurf = nodataval #apply no data value
            
    surfaceElevGrid = surfaceElevGrid.where(surfaceElevGrid < noDataSurf)  #Replace no data values with NaNs
    bedrockElevGrid = bedrockElevGrid.where(bedrockElevGrid > noDataBed)   #Replace no data values with NaNs
    return bedrockGridAlign, surfaceGridAlign

def getDriftThick(surface, bedrock, noLayers=9, plotData=True):
    driftThick = surface - bedrock
    driftThick.data = np.clip(driftThick.data, 0, np.max(driftThick.data))
    if plotData:
        driftThick.plot()
        
    layerThick = driftThick/noLayers
    
    return driftThick, layerThick
    