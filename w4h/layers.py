import datetime
import os
import pathlib

import rioxarray as rxr
import xarray as xr
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point

def merge_tables(data_df, header_df, data_cols=None, header_cols=None, on='API_NUMBER', how='inner', auto_pick_cols=False):
    if auto_pick_cols:
        header_cols = ['API_NUMBER', 'LATITUDE', 'LONGITUDE', 'BEDROCK_ELEV_FT', 'SURFACE_ELEV_FT', 'BEDROCK_DEPTH_FT', 'LAYER_THICK_FT']
        for c in header_df.columns:
            if c.startswith('ELEV_') or c.startswith('DEPTH'):
                header_cols.append(c)
            if '_PROJ' in c:
                header_cols.append(c)
        header_cols.append('geometry')
    elif header_cols is None:
        header_cols = header_df.columns
    else:
        header_cols = header_cols

    #If not specified, get all the cols
    if data_cols is None:
        data_cols = data_df.columns

    leftTable_join = data_df[data_cols]
    rightTable_join = header_df[header_cols]

    mergedTable = pd.merge(left=leftTable_join, right=rightTable_join, how=how, on=on)
    
    return mergedTable

def get_layer_depths(well_metadata, no_layers=9):
    for layer in range(0, no_layers): #For each layer
        #Make column names
        depthColName  = 'DEPTH_FT_LAYER'+str(layer+1)
        #depthMcolName = 'Depth_M_LAYER'+str(layer) 

        #Calculate depth to each layer at each well, in feet and meters
        well_metadata[depthColName]  = well_metadata['LAYER_THICK_FT'] * layer
        #headerData[depthMcolName] = headerData[depthColName] * 0.3048

    for layer in range(0, no_layers): #For each layer
        elevColName = 'ELEV_FT_LAYER'+str(layer+1)
        #elevMColName = 'ELEV_M_LAYER'+str(layer)
            
        well_metadata[elevColName]  = well_metadata['SURFACE_ELEV_FT'] - well_metadata['LAYER_THICK_FT'] * layer
        #headerData[elevMColName]  = headerData['SURFACE_ELEV_M'] - headerData['LAYER_THICK_M'] * layer
    return well_metadata

#Function to export the result of thickness of target sediments in each layer
def layer_target_thick(df, layers=9, return_all=False, export_dir=None, outfile_prefix='', use_depth=False):
    """Function to calculate thickness of target material in each layer at each well point

    Parameters
    ----------
    df : geopandas.geodataframe
        Geodataframe containing classified data, surface elevation, bedrock elevation, layer depths, geometry.
    layers : int, default=9
        Number of layers in model, by default 9
    return_all : bool, default=False
        If True, return list of original geodataframes with extra column added for target thick for each layer.
        If False, return list of geopandas.geodataframes with only essential information for each layer.
    export_dir : str or pathlib.Path, default=None
        If str or pathlib.Path, should be directory to which to export dataframes built in function.
    outfile_prefix : str, default=''
        Only used if export_dir is set. Will be used at the start of the exported filenames.

    Returns
    -------
    res_df or res : geopandas.geodataframe
        Geopandas geodataframe containing only important information needed for next stage of analysis.
    """
    df['TOP_ELEV_FT'] = df['SURFACE_ELEV_FT'] - df['TOP']
    df['BOT_ELEV_FT'] = df['SURFACE_ELEV_FT'] - df['BOTTOM']

    layerList = range(1,layers+1)
    res_list = []
    resdf_list = []
    #Generate Column names based on (looped) integers
    for layer in layerList:
        if use_depth:
            zStr = 'DEPTH'
            zColT = 'TOP'
            zColB = 'BOTTOM'
        else:
            zStr = 'ELEV'
            zColT = 'TOP_ELEV_FT'
            zColB = 'BOT_ELEV_FT'
        topCol = zStr+'_FT_LAYER'+str(layer)
        if layer != 9: #For all layers except the bottom layer....
            botCol = zStr+'_FT_LAYER'+str(layer+1) #use the layer below it to 
        else: #Otherwise, ...
            botCol = "BEDROCK_"+zStr+"_FT" #Use the (corrected) bedrock depth

        #Divide records into 4 separate categories for ease of calculation, to be joined back together later  
            #Category 1: Well interval starts above layer top, ends within model layer
            #Category 2: Well interval is entirely contained withing model layer
            #Category 3: Well interval starts within model layer, continues through bottom of model layer
            #Category 4: well interval begins and ends on either side of model layer (model layer is contained within well layer)

        #records1 = intervals that go through the top of the layer and bottom is within layer
        records1 = df.loc[(df[zColT] <= df[topCol]) & #Top of the well is above or equal to the top of the layer
                        (df[zColB] < df[botCol]) & #Bottom of the well is above (shallower depth) the top of the layer
                        (df[zColB] >= df[topCol]) & #Bottom is below the top of the layer
                        (df[zColB] >= df[zColT])].copy() #Bottom is deeper than top (should already be the case)
        records1['TARG_THICK'] = pd.DataFrame(np.round((records1.loc[:,topCol]-records1.loc[: , zColB]) * records1['TARGET'],3)).copy() #Multiply "target" (1 or 0) by length within layer            

        #records2 = entire interval is within layer
        records2 = df.loc[(df[zColT] >= df[topCol]) & #Top of the well is lower than top of the layer (deeper depth)
                    (df[zColB] <= df[botCol]) & #Bottom of the well is above bottom of the layer (shallower depth)
                    (df[zColB] >= df[zColT])].copy() #Bottom ofthe well is deeper than top (should already be the case)
        records2['TARG_THICK'] = pd.DataFrame(np.round((records2.loc[: , zColT] - records2.loc[: , zColB]) * records2['TARGET'],3)).copy()

        #records3 = intervals with top within layer and bottom of interval going through bottom of layer
        records3 = df.loc[(df[zColT] < df[botCol]) & #Top of the well is above (smaller depth) than bottom of layer
                    (df[zColB] > df[botCol]) & #Bottom of the well is below (deeper depth) than bottom of layer
                    (df[zColT] >= df[topCol]) & #Top of well is below (deeper depth) than top of layer
                    (df[zColB] >= df[zColT])].copy() #Bottom is deeper than top (should already be the case)
        records3['TARG_THICK'] = pd.DataFrame(np.round((records3.loc[: , zColT] - (records3.loc[:,botCol]))*records3['TARGET'],3)).copy()

        #records4 = interval goes through entire layer
        records4 = df.loc[(df[zColT] <= df[topCol]) & #Top of well is above (shallower depth) than top of layer
                    (df[zColB] > df[botCol]) & #Bottom of well is below (deeper depth) than bottom of layer
                    (df[zColB] >= df[zColT])].copy() #Bottom of well is below (deeper depth) than top of well
        records4['TARG_THICK'] = pd.DataFrame(np.round((records4.loc[: , topCol]-records4.loc[: , botCol]) * records4['TARGET'],3)).copy()
        
        #Put the four calculated record categories back together into single dataframe
        res = pd.concat([records1, records2, records3, records4])

        #The sign may be reversed if using depth rather than elevation
        if (res['TARG_THICK'] < 0).all():
            res['TARG_THICK'] = res['TARG_THICK'] * -1

        #Get geometrys for each unique API/well
        res_df = res.groupby(by=['API_NUMBER','LATITUDE','LONGITUDE'], as_index=False).sum(numeric_only=True)#Calculate thickness for each well interval in the layer indicated (e.g., if there are two well intervals from same well in one model layer)
        uniqInd = pd.DataFrame([v.values[0] for k, v in res.groupby('API_NUMBER').groups.items()]).loc[:,0]
        geomCol = res.loc[uniqInd, 'geometry']
        geomCol = pd.DataFrame(geomCol[~geomCol.index.duplicated(keep='first')]).reset_index()
        
        res_df['TARG_THICK_PER'] = pd.DataFrame(np.round(res_df['TARG_THICK']/res_df['LAYER_THICK_FT'],3)) #Calculate thickness as percent of total layer thickness
        res_df["LAYER"] = layer #Just to have as part of the output file, include the present layer in the file itself as a separate column
        res_df = res_df[['API_NUMBER', 'LATITUDE', 'LONGITUDE', 'LATITUDE_PROJ', 'LONGITUDE_PROJ','TOP', 'BOTTOM', 'TOP_ELEV_FT', 'BOT_ELEV_FT', 'SURFACE_ELEV_FT', topCol, botCol,'LAYER_THICK_FT','TARG_THICK', 'TARG_THICK_PER', 'LAYER']].copy() #Format dataframe for output
        res_df = gpd.GeoDataFrame(res_df, geometry=geomCol.loc[:,'geometry'])
        resdf_list.append(res_df)
        res_list.append(res)

        if isinstance(export_dir, pathlib.PurePath) or type(export_dir) is str:
            export_dir = pathlib.Path(export_dir)
            if export_dir.is_dir():
                pass
            else:
                try:
                    os.mkdir(export_dir)
                except:
                    print('Specified export directory does not exist and cannot be created. Function will continue run, but data will not be exported.')

            #Format and build export filepath
            export_dir = str(export_dir).replace('\\', '/').replace('\\'[0], '/')
            zFillDigs = len(str(len(layerList)))
            if outfile_prefix[-1]=='_':
                outfile_prefix = outfile_prefix[:-1]
            if export_dir[-1] =='/':
                export_dir = export_dir[:-1]
            nowStr = str(datetime.datetime.today().date())+'_'+str(datetime.datetime.today().hour)+'-'+str(datetime.datetime.today().minute)+'-'+str(datetime.datetime.today().second)
            outPath = export_dir+'/'+outfile_prefix+'_Lyr'+str(layer).zfill(zFillDigs)+'_'+nowStr+'.csv'

            if return_all:
                res.to_csv(outPath, sep=',', na_rep=np.nan, index_label='ID')
                res_df.to_csv(outPath, sep=',', na_rep=np.nan, index_label='ID')
            else:
                res_df.to_csv(outPath, sep=',', na_rep=np.nan, index_label='ID')

    if return_all:
        return res_list, resdf_list
    else:
        return resdf_list

#Interpolate layers to model grid
def layer_interp(points, grid, method='nearest', lin_kind='cubic', export_dir=None, targetcol='TARG_THICK_PER', lyrcol='LAYER', xcol=None, ycol=None, xcoord='x', ycoord='y', **kwargs):
    """Function to interpolate wells to model grid

    Parameters
    ----------
    points : _type_
        _description_
    grid : _type_
        _description_
    kind : _type_
        _description_
    layers : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    nnList = ['nearest', 'nearest neighbor', 'nearestneighbor','neighbor',  'nn','n']
    splineList = ['interp2d', 'interp', 'spline', 'spl', 'sp', 's']
    linList = ['linear', 'lin', 'l']
    linKindList = ['linear', 'cubic', 'quintic']
    ctList = ['clough tocher', 'clough', 'cloughtocher', 'ct', 'c']
    rbfList = ['rbf', 'radial basis', 'radial basis function', 'r', 'radial']
    #k-nearest neighbors from scikit-learn?
    #kriging? (from pykrige or maybe also from scikit-learn)
    
    if xcol is None:
        dataX = points['geometry'].x
    else:
        dataX = points[xcol]
    
    if ycol is None:
        dataY = points['geometry'].y
    else:
        dataY = points[ycol]

    layer = points[lyrcol]        
    interpVal = points[targetcol]

    X = grid[xcoord].values #Extract xcoords from grid
    Y = grid[ycoord].values #Extract ycoords from grid

    if method.lower() in nnList:
        X, Y = np.meshgrid(X, Y) #2D Grid for interpolation
        interp = interpolate.NearestNDInterpolator(list(zip(dataX, dataY)), interpVal, **kwargs)
        Z = interp(X, Y)
    elif method.lower() in linList:
        X, Y = np.meshgrid(X, Y) #2D Grid for interpolation
        interp = interpolate.LinearNDInterpolator(list(zip(dataX, dataY)), interpVal, **kwargs)
        Z = interp(X, Y)
    elif method.lower() in ctList:
        X, Y = np.meshgrid(X, Y) #2D Grid for interpolation
        interp = interpolate.CloughTocher2DInterpolator(list(zip(dataX, dataY)), interpVal, **kwargs)
        Z = interp(X, Y)    
    elif method.lower() in rbfList:
        dataXY=  np.column_stack((dataX, dataY))
        interp = interpolate.RBFInterpolator(dataXY, interpVal, **kwargs)
        print("Radial Basis Function does not work well with many well-based datasets. Consider instead specifying 'nearest', 'linear', 'spline', or 'clough tocher' for interpolation method.")
        Z = interp(np.column_stack((X.ravel(), Y.ravel()))).reshape(X.shape)
    elif method.lower() in splineList:
        if lin_kind.lower() in linKindList:
            interp = interpolate.interp2d(dataX, dataY, interpVal, kind=lin_kind, **kwargs)
        else:
            X, Y = np.meshgrid(X, Y) #2D Grid for interpolation
            print("Specified kind of interpolation (lin_kind={}) not recognized, using 'cubic'.".format(lin_kind))
            interp = interpolate.interp2d(list(zip(dataX, dataY)), interpVal, kwargs, kind='cubic', **kwargs)
        Z = interp(X, Y)
    else:
        print('Specified interpolation method not recognized, using nearest neighbor.')
        interp = interpolate.NearestNDInterpolator(list(zip(dataX, dataY)), interpVal, **kwargs)
        Z = interp(X, Y)
    

    interp_grid = xr.DataArray( #Create new datarray with new data values, but everything else the same
                data=Z,
                dims=grid.dims,
                coords=grid.coords,
                attrs=grid.attrs)
    #interp_grid=interp_grid.interpolate_na(dim=x)

    return interp_grid
