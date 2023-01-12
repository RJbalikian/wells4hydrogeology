import rioxarray as rxr
import xarray as xr
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point
import datetime

def mergeTables(leftTable, rightTable, leftCols='', rightCols='', onCol='API_NUMBER', how='inner'):
    if leftCols == '':
        leftCols = leftTable.columns
    if rightCols == '':
        rightCols = rightTable.columns

    leftTable_join = leftTable[leftCols]
    rightTable_join = rightTable[rightCols]

    mergedTable = pd.merge(left=leftTable_join, right=rightTable_join, how=how, on=onCol)
    
    return mergedTable