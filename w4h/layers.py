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
from w4h import logger_function

#Function to Merge tables
def merge_tables(data_df, header_df, data_cols=None, header_cols=None, auto_pick_cols=False, drop_duplicate_cols=True, log=False, verbose=False, **kwargs):
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
        Whether to autopick the columns from the metadata table. If True, the following column names are kept:['API_NUMBER', 'LATITUDE', 'LONGITUDE', 'BEDROCK_ELEV_FT', 'SURFACE_ELEV_FT', 'BEDROCK_DEPTH_FT', 'LAYER_THICK_FT'], by default False
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
def get_layer_depths(df_with_depths, no_layers=9, log=False):
    """Function to calculate depths and elevations of each model layer at each well based on surface elevation, bedrock elevation, and number of layers/layer thickness

    Parameters
    ----------
    df_with_depths : pandas.DataFrame
        Dataframe containing well metdata
    no_layers : int, default=9
        Number of layers. This should correlate with get_drift_thick() input parameter, if drift thickness was calculated using that function, by default 9.
    log : bool, default = False
        Whether to log inputs and outputs to log file.

    Returns
    -------
    pandas.Dataframe
        Dataframe containing new columns for depth to layers and elevation of layers.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    for layer in range(0, no_layers): #For each layer
        #Make column names
        depthColName  = 'DEPTH_FT_LAYER'+str(layer+1)
        #depthMcolName = 'Depth_M_LAYER'+str(layer) 

        #Calculate depth to each layer at each well, in feet and meters
        df_with_depths[depthColName]  = df_with_depths['LAYER_THICK_FT'] * layer
        #headerData[depthMcolName] = headerData[depthColName] * 0.3048

    for layer in range(0, no_layers): #For each layer
        elevColName = 'ELEV_FT_LAYER'+str(layer+1)
        #elevMColName = 'ELEV_M_LAYER'+str(layer)
            
        df_with_depths[elevColName]  = df_with_depths['SURFACE_ELEV_FT'] - df_with_depths['LAYER_THICK_FT'] * layer
        #headerData[elevMColName]  = headerData['SURFACE_ELEV_M'] - headerData['LAYER_THICK_M'] * layer
    return df_with_depths

#Function to export the result of thickness of target sediments in each layer
def layer_target_thick(df, layers=9, return_all=False, export_dir=None, outfile_prefix='', depth_top_col='TOP', depth_bot_col='BOTTOM', log=False):
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

    df['TOP_ELEV_FT'] = df['SURFACE_ELEV_FT'] - df[depth_top_col]
    df['BOT_ELEV_FT'] = df['SURFACE_ELEV_FT'] - df[depth_bot_col]

    layerList = range(1,layers+1)
    res_list = []
    resdf_list = []
    #Generate Column names based on (looped) integers
    for layer in layerList:
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
        records1 = df.loc[(df[zColT] >= df[topCol]) & #Top of the well is above or equal to the top of the layer
                        (df[zColB] <= df[topCol]) & # & #Bottom is below the top of the layer
                        (df[zColB] <= df[zColT])].copy() #Bottom is deeper than top (should already be the case)
        records1['TARG_THICK_FT'] = pd.DataFrame(np.round((records1.loc[:,topCol]-records1.loc[: , zColB]) * records1['TARGET'],3)).copy() #Multiply "target" (1 or 0) by length within layer            
        
        #records2 = entire interval is within layer
        records2 = df.loc[(df[zColT] <= df[topCol]) & #Top of the well is lower than top of the layer 
                    (df[zColB] >= df[botCol]) & #Bottom of the well is above bottom of the layer 
                    (df[zColB] <= df[zColT])].copy() #Bottom ofthe well is deeper than or equal to top (should already be the case)
        records2['TARG_THICK_FT'] = pd.DataFrame(np.round((records2.loc[: , zColT] - records2.loc[: , zColB]) * records2['TARGET'],3)).copy()

        #records3 = intervals with top within layer and bottom of interval going through bottom of layer
        records3 = df.loc[(df[zColT] > df[botCol]) & #Top of the well is above bottom of layer
                    (df[zColB] < df[botCol]) & #Bottom of the well is below bottom of layer
                    (df[zColT] <= df[topCol]) & #Top of well is below top of layer
                    (df[zColB] <= df[zColT])].copy() #Bottom is deeper than top (should already be the case)
        records3['TARG_THICK_FT'] = pd.DataFrame(np.round((records3.loc[: , zColT] - (records3.loc[:,botCol]))*records3['TARGET'],3)).copy()

        #records4 = interval goes through entire layer
        records4 = df.loc[(df[zColT] >= df[topCol]) & #Top of well is above top of layer
                    (df[zColB] < df[botCol]) & #Bottom of well is below bottom of layer
                    (df[zColB] <= df[zColT])].copy() #Bottom of well is below top of well
        records4['TARG_THICK_FT'] = pd.DataFrame(np.round((records4.loc[: , topCol]-records4.loc[: , botCol]) * records4['TARGET'],3)).copy()
        
        #Put the four calculated record categories back together into single dataframe
        res = pd.concat([records1, records2, records3, records4])

        #The sign may be reversed if using depth rather than elevation
        if (res['TARG_THICK_FT'] < 0).all():
            res['TARG_THICK_FT'] = res['TARG_THICK_FT'] * -1
        
        #Cannot have negative thicknesses
        res['TARG_THICK_FT'].clip(lower=0, inplace=True)
        res['LAYER_THICK_FT'].clip(lower=0, inplace=True)
        
        #Get geometrys for each unique API/well
        res_df = res.groupby(by=['API_NUMBER','LATITUDE','LONGITUDE'], as_index=False).sum(numeric_only=True)#Calculate thickness for each well interval in the layer indicated (e.g., if there are two well intervals from same well in one model layer)
        uniqInd = pd.DataFrame([v.values[0] for k, v in res.groupby('API_NUMBER').groups.items()]).loc[:,0]
        geomCol = res.loc[uniqInd, 'geometry']
        geomCol = pd.DataFrame(geomCol[~geomCol.index.duplicated(keep='first')]).reset_index()
        

        res_df['TARG_THICK_PER'] =  pd.DataFrame(np.round(res_df['TARG_THICK_FT']/res_df['LAYER_THICK_FT'],3)) #Calculate thickness as percent of total layer thickness
        #Replace np.inf and np.nans with 0
        res_df['TARG_THICK_PER'] = res_df['TARG_THICK_PER'].where(res_df['TARG_THICK_PER']!=np.inf, other=0) 
        res_df['TARG_THICK_PER'] = res_df['TARG_THICK_PER'].where(res_df['TARG_THICK_PER']!=np.nan, other=0) 

        res_df["LAYER"] = layer #Just to have as part of the output file, include the present layer in the file itself as a separate column
        res_df = res_df[['API_NUMBER', 'LATITUDE', 'LONGITUDE', 'LATITUDE_PROJ', 'LONGITUDE_PROJ','TOP', 'BOTTOM', 'TOP_ELEV_FT', 'BOT_ELEV_FT', 'SURFACE_ELEV_FT', topCol, botCol,'LAYER_THICK_FT','TARG_THICK_FT', 'TARG_THICK_PER', 'LAYER']].copy() #Format dataframe for output
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
def layer_interp(points, grid, layers=None, method='nearest', return_type='dataarray', export_dir=None, targetcol='TARG_THICK_PER', lyrcol='LAYER', xcol=None, ycol=None, xcoord='x', ycoord='y', log=False, verbose=False, **kwargs):
    """Function to interpolate results, going from points to grid data. Uses scipy.interpolate module.

    Parameters
    ----------
    points : list
        List containing pandas dataframes or geopandas geoadataframes containing the point data. Should be resDF_list output from layer_target_thick().
    grid : xr.DataArray or xr.Dataset
        Xarray dataarray with the coordinates/spatial reference of the output grid to interpolate to
    layers : int, default=None
        Number of layers for interpolation. If None, uses the length ofthe points list to determine number of layers. By default None.
    method : str, {'nearest', 'interp2d','linear', 'cloughtocher', 'radial basis function'}
        Type of interpolation to use. See scipy.interpolate N-D scattered. Values can be any of the following (also shown in "kind" column of N-D scattered section of table here: https://docs.scipy.org/doc/scipy/tutorial/interpolate.html). By default 'nearest'
    return_type : str, {'dataarray', 'dataset'}
        Type of xarray object to return, either xr.DataArray or xr.Dataset, by default 'dataarray.'
    export_dir : str or pathlib.Path, default=None
        Export directory for interpolated grids, using w4h.export_grids(). If None, does not export, by default None.
    targetcol : str, default = 'TARG_THICK_PER'
        Name of column in points containing data to be interpolated, by default 'TARG_THICK_PER'.
    lyrcol : str, default = 'Layer'
        Name of column containing layer number. Not currently used, by default 'LAYER'
    xcol : str, default = 'None'
        Name of column containing x coordinates. If None, will look for 'geometry' column, as in a geopandas.GeoDataframe. By default None
    ycol : str, default = 'None'
        Name of column containing y coordinates. If None, will look for 'geometry' column, as in a geopandas.GeoDataframe. By default None
    xcoord : str, default='x'
        Name of x coordinate in grid, used to extract x values of grid, by default 'x'
    ycoord : str, default='y'
        Name of y coordinate in grid, used to extract x values of grid, by default 'y'
    log : bool, default = True
        Whether to log inputs and outputs to log file.
    **kwargs
        Keyword arguments to be read directly into whichever scipy.interpolate function is designated by the method parameter.

    Returns
    -------
    interp_data : xr.DataArray or xr.Dataset, depending on return_type
        By default, returns an xr.DataArray object with the layers added as a new dimension called Layer. Can also specify return_type='dataset' to return an xr.Dataset with each layer as a separate variable.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    nnList = ['nearest', 'nearest neighbor', 'nearestneighbor','neighbor',  'nn','n']
    splineList = ['interp2d', 'interp2', 'interp', 'spline', 'spl', 'sp', 's']
    linList = ['linear', 'lin', 'l']
    ctList = ['clough tocher', 'clough', 'cloughtocher', 'ct', 'c']
    rbfList = ['rbf', 'radial basis', 'radial basis function', 'r', 'radial']
    #k-nearest neighbors from scikit-learn?
    #kriging? (from pykrige or maybe also from scikit-learn)
    
    X = np.round(grid[xcoord].values, 5)# #Extract xcoords from grid
    Y = np.round(grid[ycoord].values, 5)# #Extract ycoords from grid
    
    if layers is None and (type(points) is list or type(points) is dict):
        layers = len(points)

    if len(points) != layers:
        print('You have specified a different number of layers than what is iterable in the points argument. This may not work properly.')

    if verbose:
        print('Interpolating target lithology at each layer:')
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

        #layer = pts[lyrcol]        
        interpVal = pts[targetcol]

        #return points
        dataX.dropna(inplace=True)
        dataY = dataY.loc[dataX.index]
        dataY.dropna(inplace=True)
        interpVal = interpVal.loc[dataY.index]
        interpVal.dropna(inplace=True)
        
        dataX = dataX.loc[interpVal.index]
        dataY = dataY.loc[interpVal.index]
        
        dataX = dataX.reset_index(drop=True)
        dataY = dataY.reset_index(drop=True)
        interpVal = interpVal.reset_index(drop=True)
        
        if method.lower() in nnList:
            interpType = 'Nearest Neighbor'
            X, Y = np.meshgrid(X, Y, sparse=True) #2D Grid for interpolation
            dataPoints = np.array(list(zip(dataX, dataY)))
            interp = interpolate.NearestNDInterpolator(dataPoints, interpVal, **kwargs)
            Z = interp(X, Y)
        elif method.lower() in linList:
            interpType = 'Linear' 
            dataPoints = np.array(list(zip(dataX, dataY)))
            interp = interpolate.LinearNDInterpolator(dataPoints, interpVal, **kwargs)
            X, Y = np.meshgrid(X, Y, sparse=True) #2D Grid for interpolation
            Z = interp(X, Y)
        elif method.lower() in ctList:
            interpType = 'Clough-Toucher'
            X, Y = np.meshgrid(X, Y, sparse=True) #2D Grid for interpolation
            if 'tol' not in kwargs:
                kwargs['tol'] = 1e10
            interp = interpolate.CloughTocher2DInterpolator(list(zip(dataX, dataY)), interpVal, **kwargs)
            Z = interp(X, Y) 
        elif method.lower() in rbfList:
            interpType = 'Radial Basis'
            dataXY=  np.column_stack((dataX, dataY))
            interp = interpolate.RBFInterpolator(dataXY, interpVal, **kwargs)
            print("Radial Basis Function does not work well with many well-based datasets. Consider instead specifying 'nearest', 'linear', 'spline', or 'clough tocher' for interpolation method.")
            Z = interp(np.column_stack((X.ravel(), Y.ravel()))).reshape(X.shape)
        elif method.lower() in splineList:
            interpType = 'Spline Interpolation'
            Z = interpolate.bisplrep(dataX, dataY, interpVal, **kwargs)
                #interp = interpolate.interp2d(dataX, dataY, interpVal, kind=lin_kind, **kwargs)
                #Z = interp(X, Y)
        else:
            if verbose:
                print('Specified interpolation method not recognized, using nearest neighbor.')
            interpType = 'Nearest Neighbor'
            X, Y = np.meshgrid(X, Y, sparse=True) #2D Grid for interpolation
            interp = interpolate.NearestNDInterpolator(list(zip(dataX, dataY)), interpVal, **kwargs)
            Z = interp(X, Y)

        #global ZTest
        #ZTest = Z

        interp_grid = xr.DataArray( #Create new datarray with new data values, but everything else the same
                    data=Z,
                    dims=grid.dims,
                    coords=grid.coords)
        
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

        #interp_grid=interp_grid.interpolate_na(dim=x)
        zFillDigs = len(str(layers))
        daDict['Layer'+str(lyr).zfill(zFillDigs)] = interp_grid
        del interp_grid
        if verbose:
            print('\tCompleted {} interpolation for Layer {}'.format(str(interpType).lower(), str(lyr).zfill(zFillDigs)))

    dataAList = ['dataarray', 'da', 'a', 'array']
    dataSList = ['dataset', 'ds', 'set']

    if return_type.lower() in dataAList:
        interp_data = xr.concat(daDict.values(), dim='Layer')
        interp_data = interp_data.assign_coords(Layer=np.arange(1,10))
    elif return_type.lower() in dataSList:
        interp_data = xr.Dataset(daDict)
        if verbose:
            print('Done with interpolation, getting global attrs')
        common_attrs = {}
        for i, (var_name, data_array) in enumerate(interp_data.data_vars.items()):
            if i == 0:
                common_attrs = data_array.attrs
            else:
                common_attrs = {k: v for k, v in common_attrs.items() if k in data_array.attrs and data_array.attrs[k] == v}
        interp_data.attrs.update(common_attrs)
    else:
        print("{} is not a valid input for return_type. Please set return_type to either 'dataarray' or 'dataset'".format(return_type))
        return

    if verbose:
        for i, layer_data in enumerate(interp_data):
            pts = points[i]
            pts.plot(c=pts[targetcol])

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