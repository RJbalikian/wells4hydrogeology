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

#Function to export the result of thickness of target sediments in each layer
def LayerTargetThick(df, layer=1):
    #Generate Column names based on (looped) integers
    topCol = 'DEPTH_FT_LAYER'+str(layer)
    if layer != 9: #For all layers except the bottom layer....
        botCol = 'DEPTH_FT_LAYER'+str(layer+1) #use the layer below it to 
    else: #Otherwise, ...
        botCol = "BEDROCK_DEPTH_FT" #Use the (corrected) bedrock depth

    #Divide records into 4 separate categories for ease of calculation, to be joined back together later  
        #Category 1: Well interval starts above layer top, ends within model layer
        #Category 2: Well interval is entirely contained withing model layer
        #Category 3: Well interval starts within model layer, continues through bottom of model layer
        #Category 4: well interval begins and ends on either side of model layer (model layer is contained within well layer)

    #records1 = intervals that go through the top of the layer and bottom is within layer
    records1 = df.loc[(df['TOP'] <= df[topCol]) & #Top of the well is above or equal to the top of the layer
                    (df['BOTTOM'] < df[botCol]) & #Bottom of the well is above (shallower depth) the top of the layer
                    (df['BOTTOM'] >= df[topCol]) & #Bottom is below the top of the layer
                    (df['BOTTOM'] >= df['TOP'])].copy() #Bottom is deeper than top (should already be the case)
    records1['TARG_THICK'] = pd.DataFrame(np.round((records1.loc[:,topCol]-records1.loc[: , 'BOTTOM']) * records1['TARGET'],3)).copy() #Multiply "target" (1 or 0) by length within layer

    #records2 = entire interval is within layer
    records2 = df.loc[(df['TOP'] >= df[topCol]) & #Top of the well is lower than top of the layer (deeper depth)
                (df['BOTTOM'] <= df[botCol]) & #Bottom of the well is above bottom of the layer (shallower depth)
                (df['BOTTOM'] >= df['TOP'])].copy() #Bottom ofthe well is deeper than top (should already be the case)
    records2['TARG_THICK'] = pd.DataFrame(np.round((records2.loc[: , 'TOP'] - records2.loc[: , 'BOTTOM']) * records2['TARGET'],3)).copy()

    #records3 = intervals with top within layer and bottom of interval going through bottom of layer
    records3 = df.loc[(df['TOP'] < df[botCol]) & #Top of the well is above (smaller depth) than bottom of layer
                (df['BOTTOM'] > df[botCol]) & #Bottom of the well is below (deeper depth) than bottom of layer
                (df['TOP'] >= df[topCol]) & #Top of well is below (deeper depth) than top of layer
                (df['BOTTOM'] >= df['TOP'])].copy() #Bottom is deeper than top (should already be the case)
    records3['TARG_THICK'] = pd.DataFrame(np.round((records3.loc[: , 'TOP'] - (records3.loc[:,botCol]))*records3['TARGET'],3)).copy()

    #records4 = interval goes through entire layer
    records4 = df.loc[(df['TOP'] <= df[topCol]) & #Top of well is above (shallower depth) than top of layer
                (df['BOTTOM'] > df[botCol]) & #Bottom of well is below (deeper depth) than bottom of layer
                (df['BOTTOM'] >= df['TOP'])].copy() #Bottom of well is below (deeper depth) than top of well
    records4['TARG_THICK'] = pd.DataFrame(np.round((records4.loc[: , topCol]-records4.loc[: , botCol]) * records4['TARGET'],3)).copy()
    
    #Put the four calculated record categories back together into single dataframe
    res = pd.concat([records1, records2, records3, records4])

    res_df = res.groupby(by=['API_NUMBER','LATITUDE','LONGITUDE'], as_index=False).sum(numeric_only=True)#Calculate thickness for each well interval in the layer indicated (e.g., if there are two well intervals from same well in one model layer)
    
    res_df['TARG_THICK_PER'] = pd.DataFrame(np.round(res_df['TARG_THICK']/res_df['LAYER_THICK_FT'],3)) #Calculate thickness as percent of total layer thickness
    res_df["LAYER"] = layer #Just to have as part of the output file, include the present layer in the file itself as a separate column
    res_df = res_df[['API_NUMBER', 'LATITUDE', 'LONGITUDE', 'TOP', 'BOTTOM','SURFACE_ELEV_FT', topCol,botCol,'LAYER_THICK_FT','TARG_THICK', 'TARG_THICK_PER', 'LAYER']].copy() #Format dataframe for output
    
    return res, res_df