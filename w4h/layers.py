"""The Layers module contains functions for splitting data into a layered model
and for interpolating data within the layers
"""

import datetime
import inspect
import os
import pathlib

import rioxarray as rxr
import xarray as xr
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import Point
from scipy import interpolate

import w4h
from w4h import logger_function, verbose_print

#Function to Merge tables
def merge_metadata(data_df, header_df, data_cols=None, header_cols=None, auto_pick_cols=False, drop_duplicate_cols=True, log=False, verbose=False, **kwargs):
    """Function to merge tables, intended for merging metadata table with data table

    Parameters
    ----------
    data_df : pandas.DataFrame
        "Left" dataframe, intended for this purpose to be dataframe with main data, but can be anything
    header_df : pandas.DataFrame
        "Right" dataframe, intended for this purpose to be dataframe with metadata, but can be anything
    data_cols : list, optional
        List of strings of column names, for columns to be included after join from "left" table (data table). If None, all columns are kept, by default None
    header_cols : list, optional
        List of strings of columns names, for columns to be included in merged table after merge from "right" table (metadata). If None, all columns are kept, by default None
    auto_pick_cols : bool, default = False
        Whether to autopick the columns from the metadata table. If True, the following column names are kept:['API_NUMBER', 'LATITUDE', 'LONGITUDE', 'BEDROCK_ELEV', 'SURFACE_ELEV', 'BEDROCK_DEPTH', 'LAYER_THICK'], by default False
    drop_duplicate_cols : bool, optional
        If True, drops duplicate columns from the tables so that columns do not get renamed upon merge, by default True
    log : bool, default = False
        Whether to log inputs and outputs to log file.
    **kwargs
        kwargs that are passed directly to pd.merge(). By default, the 'on' and 'how' parameters are defined as on='API_NUMBER' and how='inner'

    Returns
    -------
    mergedTable : pandas.DataFrame
        Merged dataframe
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)
    if verbose:
        verbose_print(merge_metadata, locals(), exclude_params=['data_df', 'header_df'])
    if header_df is None:
        #Figure out which columns to include
        if data_cols is None:
            #If not specified, get all the cols
            data_cols = data_df.columns                      
        else:
            if header_cols is not None:
                data_cols = data_cols + header_cols
        
        data_df = data_df[data_cols]

        mergedTable = data_df
    else:   
        if auto_pick_cols:
            header_cols = ['API_NUMBER', 'LATITUDE', 'LONGITUDE', 'BEDROCK_ELEV', 'SURFACE_ELEV', 'BEDROCK_DEPTH', 'LAYER_THICK']
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

        #Defults for on and how
        if 'on' not in kwargs.keys():
            kwargs['on']='API_NUMBER'

        if 'how' not in kwargs.keys():
            kwargs['how']='inner'

        #Drop duplicate columns
        if drop_duplicate_cols:
            header_colCopy= header_cols.copy()
            remCount = 0
            for i, c in enumerate(header_colCopy):
                if c in data_cols and c != kwargs['on']:
                    if verbose:
                        print('Removing {} (duplicate columns) from data.'.format(header_cols[i-remCount]))
                    header_cols.pop(i - remCount)
                    remCount += 1

        leftTable_join = data_df[data_cols]
        rightTable_join = header_df[header_cols]

        mergedTable = pd.merge(left=leftTable_join, right=rightTable_join, **kwargs)
    return mergedTable

#Get layer depths of each layer, based on precalculated layer thickness
def get_layer_depths(df_with_depths, surface_elev_col='SURFACE_ELEV', layer_thick_col='LAYER_THICK', layers=9, log=False):
    """Function to calculate depths and elevations of each model layer at each well based on surface elevation, bedrock elevation, and number of layers/layer thickness

    Parameters
    ----------
    df_with_depths : pandas.DataFrame
        Dataframe containing well metdata
    layers : int, default=9
        Number of layers. This should correlate with get_drift_thick() input parameter, if drift thickness was calculated using that function, by default 9.
    log : bool, default = False
        Whether to log inputs and outputs to log file.

    Returns
    -------
    pandas.Dataframe
        Dataframe containing new columns for depth to layers and elevation of layers.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    for layer in range(0, layers): #For each layer
        #Make column names
        depthColName  = 'DEPTH_LAYER'+str(layer+1)
        #depthMcolName = 'Depth_M_LAYER'+str(layer) 

        #Calculate depth to each layer at each well, in feet and meters
        df_with_depths[depthColName]  = df_with_depths[layer_thick_col] * layer
        #headerData[depthMcolName] = headerData[depthColName] * 0.3048

    for layer in range(0, layers): #For each layer
        elevColName = 'ELEV_LAYER'+str(layer+1)
        #elevMColName = 'ELEV_M_LAYER'+str(layer)
            
        df_with_depths[elevColName]  = df_with_depths[surface_elev_col] - df_with_depths[layer_thick_col] * layer
        #headerData[elevMColName]  = headerData['SURFACE_ELEV_M'] - headerData['LAYER_THICK_M'] * layer
        
    return df_with_depths

#Function to export the result of thickness of target sediments in each layer
def layer_target_thick(df, layers=9, return_all=False, export_dir=None, outfile_prefix=None, depth_top_col='TOP', depth_bot_col='BOTTOM', log=False):
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
    outfile_prefix : str, default=None
        Only used if export_dir is set. Will be used at the start of the exported filenames
    depth_top_col : str, default='TOP'
        Name of column containing data for depth to top of described well intervals
    depth_bot_col : str, default='BOTTOM'
        Name of column containing data for depth to bottom of described well intervals
    log : bool, default = True
        Whether to log inputs and outputs to log file.
    
    Returns
    -------
    res_df or res : geopandas.geodataframe
        Geopandas geodataframe containing only important information needed for next stage of analysis.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)
    df['TOP_ELEV'] = df['SURFACE_ELEV'] - df[depth_top_col]
    df['BOT_ELEV'] = df['SURFACE_ELEV'] - df[depth_bot_col]

    layerList = range(1,layers+1)
    res_list = []
    resdf_list = []

    #Generate Column names based on (looped) integers
    for layer in layerList:
        zStr = 'ELEV'
        wellIntTop_col = 'TOP_ELEV'
        wellIntBot_col = 'BOT_ELEV'
        modelLyrTop_Col = zStr+'_LAYER'+str(layer)
        if layer != 9: #For all layers except the bottom layer....
            modelLyrBot_col = zStr+'_LAYER'+str(layer+1) #use the layer below it to 
        else: #Otherwise, ...
            modelLyrBot_col = "BEDROCK_"+zStr #Use the (corrected) bedrock depth

        #Divide records into 4 categories for ease of calculation, to be joined back together later  
            #Category 1: Well interval starts above layer top, ends within model layer
            #Category 2: Well interval is entirely contained withing model layer
            #Category 3: Well interval starts within model layer, continues through bottom of model layer
            #Category 4: well interval begins and ends on either side of model layer (model layer is contained within well layer)
        df['API_NUMBER'] = df['API_NUMBER'].astype(int).astype(str)
        #records1 = intervals that go through the top of the layer and bottom is within layer
        records1 = df.loc[(df[wellIntTop_col] > df[modelLyrTop_Col]) & #Top of the well interval is above or equal to the top of the layer
                        (df[wellIntBot_col] <= df[modelLyrTop_Col]) & # & #Bottom is below the top of the layer
                        (df[wellIntBot_col] <= df[wellIntTop_col])].copy() #Bottom is deeper than top (should already be the case)
        records1['TARG_THICK'] = pd.DataFrame(np.round((records1.loc[:,modelLyrTop_Col]-records1.loc[: , wellIntBot_col]) * records1['TARGET'],3)).copy() #Multiply "target" (1 or 0) by length within layer            
        
        #records2 = entire interval is within layer
        records2 = df.loc[(df[wellIntTop_col] <= df[modelLyrTop_Col]) & #Top of the well is lower than top of the layer 
                    (df[wellIntBot_col] >= df[modelLyrBot_col]) & #Bottom of the well is above bottom of the layer 
                    (df[wellIntBot_col] <= df[wellIntTop_col])].copy() #Bottom ofthe well is deeper than or equal to top (should already be the case)
        records2['TARG_THICK'] = pd.DataFrame(np.round((records2.loc[: , wellIntTop_col] - records2.loc[: , wellIntBot_col]) * records2['TARGET'],3)).copy()

        #records3 = intervals with top within layer and bottom of interval going through bottom of layer
        records3 = df.loc[(df[wellIntTop_col] >= df[modelLyrBot_col]) & #Top of the well is above bottom of layer
                    (df[wellIntTop_col] <= df[modelLyrTop_Col]) & #Top of well is below top of layer
                    (df[wellIntBot_col] < df[modelLyrBot_col]) & #Bottom of the well is below bottom of layer
                    (df[wellIntBot_col] <= df[wellIntTop_col])].copy() #Bottom is deeper than top (should already be the case)
        records3['TARG_THICK'] = pd.DataFrame(np.round((records3.loc[: , wellIntTop_col] - (records3.loc[:,modelLyrBot_col]))*records3['TARGET'],3)).copy()

        #records4 = interval goes through entire layer
        records4 = df.loc[(df[wellIntTop_col] > df[modelLyrTop_Col]) & #Top of well is above top of layer
                    (df[wellIntBot_col] < df[modelLyrBot_col]) & #Bottom of well is below bottom of layer
                    (df[wellIntBot_col] <= df[wellIntTop_col])].copy() #Bottom of well is below top of well
        records4['TARG_THICK'] = pd.DataFrame(np.round((records4.loc[: , modelLyrTop_Col]-records4.loc[: , modelLyrBot_col]) * records4['TARGET'],3)).copy()

        # Truth check
        inputRecordList = [records1, records3, records4]
        truthCheckedRecords = []
        for i, recDF in enumerate(inputRecordList):
            if recDF.shape[0]>1:
                recDF = recDF.loc[recDF['TARG_THICK'].idxmax()]
            
            truthCheckedRecords.append(recDF)

        # Truth check records category 2, which may have multiples
        rec2Sorted = records2.sort_values(by='TARG_THICK', ascending=False)
        subS = ['API_NUMBER', 'TOP', "BOTTOM"]
        records2 = rec2Sorted.drop_duplicates(subset=subS, 
                                            keep='first')

        truthCheckedRecords.insert(1, records2)
        
        #Put the four calculated record categories back together into single dataframe
        res = pd.concat(truthCheckedRecords)

        #The sign may be reversed if using depth rather than elevation
        if (res['TARG_THICK'] < 0).all():
            res['TARG_THICK'] = res['TARG_THICK'] * -1
        
        #Cannot have negative thicknesses
        res['TARG_THICK'] = res['TARG_THICK'].clip(lower=0)
        res['LAYER_THICK'] = res['LAYER_THICK'].clip(lower=0)

        #Get geometries for each unique API/well
        res.reset_index(drop=True, inplace=True)
        #Calculate thickness for each well interval in the layer indicated (e.g., if there are two well intervals from same well in one model layer)
        res_df = res.groupby(by=['API_NUMBER','LATITUDE','LONGITUDE'], as_index=False)['TARG_THICK'].sum()

        mergeRes = res.drop(['LATITUDE', "LONGITUDE", "TARG_THICK"], axis=1)
        res_df = pd.merge(left=res_df, right=mergeRes, how='left', on='API_NUMBER')
        
        testDF = pd.DataFrame()
        testDF['TARG_THICK'] = res_df['TARG_THICK'].clip(upper=res_df['LAYER_THICK'])

        
        uniqInd = pd.DataFrame([v.values[0] for k, v in res.groupby(by=['API_NUMBER','LATITUDE','LONGITUDE']).groups.items()]).loc[:,0]
        geomCol = res.loc[uniqInd, 'geometry']
        geomCol = pd.DataFrame(geomCol[~geomCol.index.duplicated(keep='first')]).reset_index()
        
        #Calculate thickness as percent of total layer thickness
        res_df['TARG_THICK_PER'] =  pd.DataFrame(np.round(res_df['TARG_THICK']/res_df['LAYER_THICK'],3)).clip(upper=1)
        #Replace np.inf and np.nans with 0
        res_df['TARG_THICK_PER'] = res_df['TARG_THICK_PER'].where(res_df['TARG_THICK_PER']!=np.inf, other=0)
        res_df['TARG_THICK_PER'] = res_df['TARG_THICK_PER'].where(res_df['TARG_THICK_PER']!=np.nan, other=0)

        res_df["LAYER"] = layer #Just to have as part of the output file, include the present layer in the file itself as a separate column
        res_df = res_df[['API_NUMBER', 'LATITUDE', 'LONGITUDE', 'LATITUDE_PROJ', 'LONGITUDE_PROJ','TOP', 'BOTTOM', 'TOP_ELEV', 'BOT_ELEV', 'SURFACE_ELEV', modelLyrTop_Col, modelLyrBot_col,'LAYER_THICK','TARG_THICK', 'TARG_THICK_PER', 'LAYER']].copy() #Format dataframe for output
        res_df = gpd.GeoDataFrame(res_df, geometry=geomCol.loc[:,'geometry'])

        resdf_list.append(res_df)
        res_list.append(res)

        if isinstance(export_dir, (pathlib.PurePath, str)):
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
            if str(outfile_prefix).endswith('_'):
                outfile_prefix = outfile_prefix[:-1]
            elif outfile_prefix is None:
                outfile_prefix = ''
            

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

#Interpolate layers to model model_grid
def layer_interp(points, model_grid, layers=None, interp_kind='nearest', surface_grid=None, bedrock_grid=None, layer_thick_grid=None, drift_thick_grid=None, return_type='dataset', export_dir=None, target_col='TARG_THICK_PER', layer_col='LAYER', xcol=None, ycol=None, xcoord='x', ycoord='y', log=False, verbose=False, **kwargs):
    """Function to interpolate results, going from points to model_grid data. Uses scipy.interpolate module.

    Parameters
    ----------
    points : list
        List containing pandas dataframes or geopandas geoadataframes containing the point data. Should be resDF_list output from layer_target_thick().
    model_grid : xr.DataArray or xr.Dataset
        Xarray DataArray or DataSet with the coordinates/spatial reference of the output model_grid to interpolate to
    layers : int, default=None
        Number of layers for interpolation. If None, uses the length ofthe points list to determine number of layers. By default None.
    interp_kind : str, {'nearest', 'interp2d','linear', 'cloughtocher', 'radial basis function'}
        Type of interpolation to use. See scipy.interpolate N-D scattered. Values can be any of the following (also shown in "kind" column of N-D scattered section of table here: https://docs.scipy.org/doc/scipy/tutorial/interpolate.html). By default 'nearest'
    return_type : str, {'dataset', 'dataarray'}
        Type of xarray object to return, either xr.DataArray or xr.Dataset, by default 'dataset.'
    export_dir : str or pathlib.Path, default=None
        Export directory for interpolated grids, using w4h.export_grids(). If None, does not export, by default None.
    target_col : str, default = 'TARG_THICK_PER'
        Name of column in points containing data to be interpolated, by default 'TARG_THICK_PER'.
    layer_col : str, default = 'Layer'
        Name of column containing layer number. Not currently used, by default 'LAYER'
    xcol : str, default = 'None'
        Name of column containing x coordinates. If None, will look for 'geometry' column, as in a geopandas.GeoDataframe. By default None
    ycol : str, default = 'None'
        Name of column containing y coordinates. If None, will look for 'geometry' column, as in a geopandas.GeoDataframe. By default None
    xcoord : str, default='x'
        Name of x coordinate in model_grid, used to extract x values of model_grid, by default 'x'
    ycoord : str, default='y'
        Name of y coordinate in model_grid, used to extract x values of model_grid, by default 'y'
    log : bool, default = True
        Whether to log inputs and outputs to log file.
    **kwargs
        Keyword arguments to be read directly into whichever scipy.interpolate function is designated by the interp_kind parameter.

    Returns
    -------
    interp_data : xr.DataArray or xr.Dataset, depending on return_type
        By default, returns an xr.DataArray object with the layers added as a new dimension called Layer. Can also specify return_type='dataset' to return an xr.Dataset with each layer as a separate variable.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if verbose:
        verbose_print(layer_interp, locals(), exclude_params=['points', 'model_grid'])

    nnList = ['nearest', 'nearest neighbor', 'nearestneighbor','neighbor',  'nn','n']
    splineList = ['interp2d', 'interp2', 'interp', 'spline', 'spl', 'sp', 's']
    linList = ['linear', 'lin', 'l']
    ctList = ['clough tocher', 'clough', 'cloughtocher', 'ct', 'c']
    rbfList = ['rbf', 'radial basis', 'radial basis function', 'r', 'radial']
    #k-nearest neighbors from scikit-learn?
    #kriging? (from pykrige or maybe also from scikit-learn)
    
    X = np.round(model_grid[xcoord].values, 5)  # Extract xcoords from model_grid
    Y = np.round(model_grid[ycoord].values, 5)  # Extract ycoords from model_grid
    
    if layers is None and (type(points) is list or type(points) is dict):
        layers = len(points)

    if len(points) != layers:
        print('You have specified a different number of layers than what is iterable in the points argument. This may not work properly.')

    if verbose:
        print('\tInterpolating target lithology at each layer:')
    
    daDict = {}
    for lyr in range(1, layers+1):
        
        if type(points) is list or type(points) is dict:
            pts = points[lyr-1]
        else:
            pts = points

        if xcol is None:
            if 'geometry' in pts.columns:
                dataX = pts['geometry'].x
            else:
                print('xcol not specified and geometry column not detected (points is not/does not contain geopandas.GeoDataFrame)' )
                return
        else:
            dataX = pts[xcol]
        
        if ycol is None:
            if 'geometry' in pts.columns:
                dataY = pts['geometry'].y
            else:
                print('ycol not specified and geometry column not detected (points is not/does not contain geopandas.GeoDataFrame)' )
                return
        else:
            dataY = pts[ycol]

        interpVal = pts[target_col]

        # Return points
        dataX = dataX.dropna()
        dataY = dataY.loc[dataX.index]
        dataY = dataY.dropna()
        interpVal = interpVal.loc[dataY.index]
        interpVal = interpVal.dropna()
        
        dataX = dataX.loc[interpVal.index]
        dataY = dataY.loc[interpVal.index]
        
        dataX = dataX.reset_index(drop=True)
        dataY = dataY.reset_index(drop=True)
        interpVal = interpVal.reset_index(drop=True)
        
        if interp_kind.lower() in nnList:
            interpType = 'Nearest Neighbor'
            X, Y = np.meshgrid(X, Y, sparse=True) #2D Grid for interpolation
            dataPoints = np.array(list(zip(dataX, dataY)))
            interp = interpolate.NearestNDInterpolator(dataPoints, interpVal, **kwargs)
            Z = interp(X, Y)
        elif interp_kind.lower() in linList:
            interpType = 'Linear' 
            dataPoints = np.array(list(zip(dataX, dataY)))
            interp = interpolate.LinearNDInterpolator(dataPoints, interpVal, **kwargs)
            X, Y = np.meshgrid(X, Y, sparse=True) #2D Grid for interpolation
            Z = interp(X, Y)
        elif interp_kind.lower() in ctList:
            interpType = 'Clough-Toucher'
            X, Y = np.meshgrid(X, Y, sparse=True) #2D Grid for interpolation
            if 'tol' not in kwargs:
                kwargs['tol'] = 1e10
            interp = interpolate.CloughTocher2DInterpolator(list(zip(dataX, dataY)), interpVal, **kwargs)
            Z = interp(X, Y) 
        elif interp_kind.lower() in rbfList:
            interpType = 'Radial Basis'
            dataXY=  np.column_stack((dataX, dataY))
            interp = interpolate.RBFInterpolator(dataXY, interpVal, **kwargs)
            print("Radial Basis Function does not work well with many well-based datasets. Consider instead specifying 'nearest', 'linear', 'spline', or 'clough tocher' for interpolation interp_kind.")
            Z = interp(np.column_stack((X.ravel(), Y.ravel()))).reshape(X.shape)
        elif interp_kind.lower() in splineList:
            interpType = 'Spline Interpolation'
            Z = interpolate.bisplrep(dataX, dataY, interpVal, **kwargs)
                #interp = interpolate.interp2d(dataX, dataY, interpVal, kind=lin_kind, **kwargs)
                #Z = interp(X, Y)
        else:
            if verbose:
                print('Specified interpolation interp_kind not recognized, using nearest neighbor.')
            interpType = 'Nearest Neighbor'
            X, Y = np.meshgrid(X, Y, sparse=True) #2D Grid for interpolation
            interp = interpolate.NearestNDInterpolator(list(zip(dataX, dataY)), interpVal, **kwargs)
            Z = interp(X, Y)

        interp_grid = xr.DataArray( #Create new datarray with new data values, but everything else the same
                    data=Z,
                    dims=model_grid.dims,
                    coords=model_grid.coords)
        
        if 'band' in interp_grid.coords:
            interp_grid = interp_grid.drop_vars('band')
        interp_grid = interp_grid.clip(min=0, max=1, keep_attrs=True)

        interp_grid = interp_grid.expand_dims(dim='Layer')
        interp_grid = interp_grid.assign_coords(Layer=[lyr])

        del Z
        del dataX
        del dataY
        del interpVal
        del interp

        zFillDigs = len(str(layers))
        daDict['Layer'+str(lyr).zfill(zFillDigs)] = interp_grid
        del interp_grid
        if verbose:
            print('\t\tCompleted {} interpolation for Layer {}'.format(str(interpType).lower(), str(lyr).zfill(zFillDigs)))



    dataArrayList = ['dataarray', 'da', 'a', 'array']
    dataSetList = ['dataset', 'ds', 'set']

    interp_data = xr.concat(daDict.values(), dim='Layer')
    interp_data = interp_data.assign_coords(Layer=np.arange(1, layers+1))
    
    if return_type.lower() in dataArrayList:
        pass        
    else:
        if return_type.lower() not in dataSetList and verbose:
            print(f"{return_type} is not a valid input for return_type. Please set return_type to either 'dataarray' or 'dataset'")
            print("Using dataset by default")            
        interp_data = xr.Dataset(data_vars={'Model_Layers': interp_data})
        if verbose:
            print('Done with interpolation, getting additional layers and attributes')

        # Get common attributes from all layers to use as "global" attributes
        common_attrs = {}
        for i, (var_name, data_array) in enumerate(interp_data.data_vars.items()):
            if i == 0:
                common_attrs = data_array.attrs
            else:
                common_attrs = {k: v for k, v in common_attrs.items() if k in data_array.attrs and data_array.attrs[k] == v}
        interp_data.attrs.update(common_attrs)

    if verbose:
        for i, layer_data in enumerate(interp_data):
            pts = points[i]
            pts.plot(c=pts[target_col])

    if export_dir is None:
        pass
    else:
        w4h.export_grids(grid_data=interp_data, out_path=export_dir, file_id='',filetype='tif', variable_sep=True, date_stamp=True)
        print('Exported to {}'.format(export_dir))

    return interp_data

#Optional, combine dataset
def combine_dataset(layer_dataset, surface_elev, bedrock_elev, layer_thick, log=False):
    """Function to combine xarray datasets or datarrays into a single xr.Dataset. Useful to add surface, bedrock, layer thick, and layer datasets all into one variable, for pickling, for example.

    Parameters
    ----------
    layer_dataset : xr.DataArray 
        DataArray contining all the interpolated layer information.
    surface_elev : xr.DataArray
        DataArray containing surface elevation data
    bedrock_elev : xr.DataArray
        DataArray containing bedrock elevation data
    layer_thick : xr.DataArray
        DataArray containing the layer thickness at each point in the model grid
    log : bool, default = False
        Whether to log inputs and outputs to log file.

    Returns
    -------
    xr.Dataset
        Dataset with all input arrays set to different variables within the dataset.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)
    
    daDict = {}
    daDict['Layers'] = layer_dataset
    daDict['Surface_Elev'] = surface_elev
    daDict['Bedrock_Elev'] = bedrock_elev
    daDict['Layer_Thickness'] = layer_thick

    combined_dataset = xr.Dataset(daDict)

    return combined_dataset
