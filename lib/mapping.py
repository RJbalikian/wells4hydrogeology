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

def readWCS():

    wcs_url = r'https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_DEM_30M/ImageServer/WCSServer?request=GetCapabilities&service=WCS'

    # Create coverage object
    my_wcs = WebCoverageService(wcs_url,
                                 version='2.0.1')

    # Get list of coverages
    print(my_wcs.contents.keys())

    # Get geo-bounding boxes and native CRS
    #my_wcs.contents['AverageChlorophyllScaled'].boundingboxes

    # Get axis labels
    #my_wcs.contents['AverageChlorophyllScaled'].grid.axislabels

    # Get dimension
    #my_wcs.contents['AverageChlorophyllScaled'].grid.dimension

    # Get grid lower and upper bounds
    #my_wcs.contents['AverageChlorophyllScaled'].grid.lowlimits

    #my_wcs.contents['AverageChlorophyllScaled'].grid.highlimits

    #my_wcs.contents['AverageChlorophyllScaled'].grid.offsetvectors

    # For coverage with time axis get the date time values
    #my_wcs.contents['AverageChlorophyllScaled'].timepositions
    return