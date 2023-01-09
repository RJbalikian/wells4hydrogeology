import rioxarray
import geopandas as gpd
import pandas as pd
import numpy as np
from owslib.wcs import WebCoverageService

def rastertoPoints_extract(raster, points):
    
    return points

def addElevtoHeader(xyz, header):
    headerXYZData = header.merge(xyz, how='left', on='API_NUMBER')
    headerXYZData.drop(['LATITUDE_x', 'LONGITUDE_x'], axis=1, inplace=True)
    headerXYZData.rename({'LATITUDE_y':'LATITUDE', 'LONGITUDE_y':'LONGITUDE'}, axis=1, inplace=True)
    return headerXYZData


def readWCS(wcs_url):

    #wcs_url = r"https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_DEM_30M/ImageServer/WCSServer?request=GetCapabilities&service=WCS"

    # Create coverage object
    my_wcs = WebCoverageService(wcs_url,
                                version='1.1.2')

    # Get list of coverages
    print(my_wcs.contents.keys())

    #minLat = 39.991998
    #maxLat = 40.248993
    #minlon = -88.493952
    #maxLon = -88.079875
    #(testWCS.getCoverage())
    return my_wcs