"""The Layers module contains functions for splitting data into a layered model
and for interpolating data within the layers
"""

import datetime
import inspect
import numbers
import os
import pathlib

import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely import Point, Polygon

import w4h
from w4h import logger_function, verbose_print


# Function to Merge tables
def merge_metadata(data_df, header_df,
                   well_id_col='API_NUMBER', data_cols=None,
                   header_cols=None, auto_pick_cols=False,
                   drop_duplicate_cols=True,
                   log=False, verbose=False, **kwargs):

    """Function to merge tables, intended for merging metadata and data tables

    Parameters
    ----------
    data_df : pandas.DataFrame
        "Left" dataframe, intended for this purpose to be
        dataframe with main data, but can be anything
    header_df : pandas.DataFrame
        "Right" dataframe, intended for this purpose to be
        dataframe with metadata, but can be anything
    data_cols : list, optional
        List of strings of column names, for columns
        to be included after join from "left" table (data table).
        If None, all columns are kept, by default None
    header_cols : list, optional
        List of strings of columns names, for columns to be included
        in merged table after merge from "right" table (metadata).
        If None, all columns are kept, by default None
    auto_pick_cols : bool, default = False
        Whether to autopick the columns from the metadata table.
        If True, the following column names are kept:
        [`well_id_col`, 'LATITUDE', 'LONGITUDE', 'BEDROCK_ELEV',
        'SURFACE_ELEV', 'BEDROCK_DEPTH', 'LAYER_THICK'], by default False.
    drop_duplicate_cols : bool, optional
        If True, drops duplicate columns from the tables so that
        columns do not get renamed upon merge, by default True.
    log : bool, default = False
        Whether to log inputs and outputs to log file.
    **kwargs
        kwargs that are passed directly to pd.merge().
        By default, the 'on' and 'how' parameters are defined as
        on=`well_id_col` and how='inner'.

    Returns
    -------
    mergedTable : pandas.DataFrame
        Merged dataframe
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)
    if verbose:
        verbose_print(merge_metadata, locals(), exclude_params=['data_df',
                                                                'header_df'])
    if header_df is None:
        # Figure out which columns to include
        if data_cols is None:
            # If not specified, get all the cols
            data_cols = data_df.columns
        else:
            if header_cols is not None:
                data_cols = data_cols + header_cols

        data_df = data_df[data_cols]

        mergedTable = data_df
    else:
        if auto_pick_cols:
            header_cols = [well_id_col, 'LATITUDE', 'LONGITUDE',
                           'BEDROCK_ELEV', 'SURFACE_ELEV', 'BEDROCK_DEPTH',
                           'LAYER_THICK']

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

        # If not specified, get all the cols
        if data_cols is None:
            data_cols = data_df.columns

        # Defults for on and how
        if 'on' not in kwargs.keys():
            kwargs['on'] = well_id_col

        if 'how' not in kwargs.keys():
            kwargs['how'] = 'inner'

        # Drop duplicate columns
        if drop_duplicate_cols:
            header_colCopy = header_cols.copy()
            remCount = 0
            for i, c in enumerate(header_colCopy):
                if c in data_cols and c != kwargs['on']:
                    if verbose:
                        print('Removing {} (duplicate columns) from data.'.format(header_cols[i-remCount]))
                    header_cols.pop(i - remCount)
                    remCount += 1

        leftTable_join = data_df[data_cols]
        rightTable_join = header_df[header_cols]

        mergedTable = pd.merge(left=leftTable_join, right=rightTable_join,
                               **kwargs)
    return mergedTable


# Get layer depths of each layer, based on precalculated layer thickness
def get_layer_depths(df_with_depths, surface_elev_col='SURFACE_ELEV',
                     layer_thick_col='LAYER_THICK', layers=9, log=False):

    """Function to calculate depths and elevations of each model layer
    at each well based on surface elevation, bedrock elevation,
    and number of layers/layer thickness

    Parameters
    ----------
    df_with_depths : pandas.DataFrame
        DataFrame containing well metdata
    layers : int, default=9
        Number of layers. This should correlate with
        get_drift_thick() input parameter, if drift thickness was calculated
        using that function, by default 9.
    log : bool, default = False
        Whether to log inputs and outputs to log file.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing new columns for depth to/elevation of layers.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    for layer in range(0, layers):  # For each layer
        # Make column names
        depthColName = 'DEPTH_LAYER'+str(layer+1)

        # Calculate depth to each layer at each well, in feet and meters
        df_with_depths[depthColName] = df_with_depths[layer_thick_col] * layer

    for layer in range(0, layers):  # For each layer
        elevColName = 'ELEV_LAYER'+str(layer+1)
        # elevMColName = 'ELEV_M_LAYER'+str(layer)

        df_with_depths[elevColName] = df_with_depths[surface_elev_col] - df_with_depths[layer_thick_col] * layer

    return df_with_depths


# Function to export the result of thickness of target sediments in each layer
def layer_target_thick(gdf, layers=9, well_id_col='API_NUMBER',
                       return_all=False, export_dir=None, outfile_prefix=None,
                       depth_top_col='TOP', depth_bot_col='BOTTOM', log=False,
                       **kwargs):
    
    """Function to calculate thickness of target material
       in each layer at each well point.
       This function loops through each model layer and divides up
       the well intervals from the input GeoDataFrame into 4 categories:
       * 1) Intervals that pierce top of model layer but end within it
       * 2) Intervals contained entirely within model layer
       * 3) Intervals that begin within the model layer but pierce the bottom
       * 4) Intervals that begin above and end below the model layer

       For each category, the amount/thickness of the target material
       within each layer is calculated. These records are then
       "truth-checked" to ensure there are not duplicates and combined.
       The percent thickness (target thickness/layer thickness) is calculated,
       before data is returned and/or exported.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        Geodataframe containing classified data, surface elevation,
        bedrock elevation, layer depths, geometry.
    layers : int, default=9
        Number of layers in model, by default 9
    well_id_col : str, default="API_NUMBER"
        The name of the column that is used for uniquely identifying each well
    return_all : bool, default=False
        If True, return list of original GeoDataFrames with extra column
        added for target thick for each layer.
        If False, return list of geopandas.GeoDataFrames with only
        essential information for each layer.
    export_dir : str or pathlib.Path, default=None
        If str or pathlib.Path, should be directory to which
        to export dataframes built in function.
    outfile_prefix : str, default=None
        Only used if export_dir is set.
        Will be used at the start of the exported filenames.
    depth_top_col : str, default='TOP'
        Name of column containing data for depth
        to top of described well intervals.
    depth_bot_col : str, default='BOTTOM'
        Name of column containing data for depth
        to bottom of described well intervals.
    log : bool, default = True
        Whether to log inputs and outputs to log file.

    Returns
    -------
    res_df and/or res : list
        A list of Geopandas GeoDataFrames containing only important information
        needed for next stage of analysis. If `return_all=True`, the input data
        with the actual descriptions will be returned
        as a separate list of GeoDataFrames.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)
    gdf['TOP_ELEV'] = gdf['SURFACE_ELEV'] - gdf[depth_top_col]
    gdf['BOT_ELEV'] = gdf['SURFACE_ELEV'] - gdf[depth_bot_col]

    layerList = range(1, layers+1)
    res_list = []
    resdf_list = []

    # Generate Column names based on (looped) integers
    for layer in layerList:
        zStr = 'ELEV'
        wellIntTop_col = 'TOP_ELEV'
        wellIntBot_col = 'BOT_ELEV'
        modelLyrTop_Col = zStr+'_LAYER'+str(layer)
        if layer != 9:  # For all layers except the bottom layer....
            modelLyrBot_col = zStr+'_LAYER'+str(layer+1)  # use the layer below
        else:  # Otherwise, ...
            modelLyrBot_col = "BEDROCK_"+zStr  # Use (corrected) brock depth

        # Divide records into 4 categories for ease of calculation:
            # 1: Well interval starts above layer top, ends within model layer
            # 2: Well interval is entirely contained withing model layer
            # 3: Well interval starts in layer, continues through bottom
            # 4: well interval begins/ends above/below (respective) model layer
        if isinstance(gdf[well_id_col], numbers.Number) or str(gdf[well_id_col]).isnumeric():
            gdf[well_id_col] = gdf[well_id_col].astype(int).astype(str)

        # records1 = intervals that go through the top of the layer and bottom is within layer
        records1 = gdf.loc[(gdf[wellIntTop_col] > gdf[modelLyrTop_Col]) &  # Top of the well interval is above or equal to the top of the layer
                           (gdf[wellIntBot_col] <= gdf[modelLyrTop_Col]) &  # Bottom is below the top of the layer
                           (gdf[wellIntBot_col] <= gdf[wellIntTop_col])].copy()  # Bottom is deeper than top (should already be the case)
        records1['TARG_THICK'] = pd.DataFrame(np.round((records1.loc[:, modelLyrTop_Col] - records1.loc[:, wellIntBot_col]) * records1['TARGET'], 3)).copy()

        # records2 = entire interval is within layer
        records2 = gdf.loc[(gdf[wellIntTop_col] <= gdf[modelLyrTop_Col]) &  # Top of the well is lower than top of the layer
                           (gdf[wellIntBot_col] >= gdf[modelLyrBot_col]) &  # Bottom of the well is above bottom of the layer
                           (gdf[wellIntBot_col] <= gdf[wellIntTop_col])].copy()  # Bottom ofthe well is deeper than or equal to top (should already be the case)
        records2['TARG_THICK'] = pd.DataFrame(np.round((records2.loc[:, wellIntTop_col] - records2.loc[:, wellIntBot_col]) * records2['TARGET'], 3)).copy()

        # records3 = intervals with top within layer and bottom of interval going through bottom of layer
        records3 = gdf.loc[(gdf[wellIntTop_col] >= gdf[modelLyrBot_col]) &  # Top of the well is above bottom of layer
                           (gdf[wellIntTop_col] <= gdf[modelLyrTop_Col]) &  # Top of well is below top of layer
                           (gdf[wellIntBot_col] < gdf[modelLyrBot_col]) &  # Bottom of the well is below bottom of layer
                           (gdf[wellIntBot_col] <= gdf[wellIntTop_col])].copy()  # Bottom is deeper than top (should already be the case)
        records3['TARG_THICK'] = pd.DataFrame(np.round((records3.loc[:, wellIntTop_col] - (records3.loc[:, modelLyrBot_col])) * records3['TARGET'], 3)).copy()

        # records4 = interval goes through entire layer
        records4 = gdf.loc[(gdf[wellIntTop_col] > gdf[modelLyrTop_Col]) &  # Top of well is above top of layer
                           (gdf[wellIntBot_col] < gdf[modelLyrBot_col]) &  # Bottom of well is below bottom of layer
                           (gdf[wellIntBot_col] <= gdf[wellIntTop_col])].copy()  # Bottom of well is below top of well
        records4['TARG_THICK'] = pd.DataFrame(np.round((records4.loc[:, modelLyrTop_Col] - records4.loc[:, modelLyrBot_col]) * records4['TARGET'], 3)).copy()

        # Truth check: ensure only well data is not being duplicated
        # Do the "partial" records first (the records that don't go through the entire layer)
        inputRecordList = [records1, records3, records4]
        truthCheckedRecords = []  # initialize output list for merging later
        subSetCols = [well_id_col, 'TOP', "BOTTOM"]
        
        # Go through each set of records, and remove duplicate well intervals
        for recDF in inputRecordList:
            # First sort by target thickness (if duplicates, keep largest)
            recDFSorted = recDF.sort_values(by='TARG_THICK', ascending=False)

            # Drop duplicates and keep largest
            recDFSorted = recDFSorted.drop_duplicates(subset=subSetCols,
                                                      keep='first')
            
            # Re-sort by the well id column for consistency's sake
            recDF = recDFSorted.sort_values(by=well_id_col)

            # Append this category to list for later merging
            truthCheckedRecords.append(recDF)

        # Truth check records category 2, which may have multiples
        # First sort by target thickness (if duplicates, keep largest)
        rec2Sorted = records2.sort_values(by='TARG_THICK', ascending=False)
        
        # Drop duplicates and keep largest
        records2 = rec2Sorted.drop_duplicates(subset=subSetCols,
                                              keep='first')
        
        # Re-sort by the well id column for consistency's sake
        records2 = records2.sort_values(by=well_id_col)

        # Insert this category to list for merging
        truthCheckedRecords.insert(1, records2)

        # Merge four record categories back into single (geo)dataframe
        res = gpd.GeoDataFrame(pd.concat(truthCheckedRecords),
                               geometry='geometry', crs=gdf.crs).sort_values(by=well_id_col)

        # The sign may be reversed if using depth rather than elevation
        if (res['TARG_THICK'] < 0).all():
            res['TARG_THICK'] = res['TARG_THICK'] * -1
        
        # Cannot have negative thicknesses
        res['TARG_THICK'] = res['TARG_THICK'].clip(lower=0)
        res['LAYER_THICK'] = res['LAYER_THICK'].clip(lower=0)

        res.reset_index(drop=True, inplace=True)
        
        # Get geometries for each unique API/well
        # Calculate thickness for each well interval in the layer indicated
        #   (e.g., if there are two well intervals from same well in one model layer)

        # Group by well id column and location, and sum TARG_THICK for each entry
        res_df = res.groupby(by=[well_id_col, 'LATITUDE', 'LONGITUDE'],
                             as_index=False)['TARG_THICK'].sum()

        # Remove the columns so they are not repeated after merge
        mergeRes = res.drop(['LATITUDE', "LONGITUDE", "TARG_THICK"], axis=1)
        
        # Merge the grouped/summed dataframe with the main dataframe
        res_df = pd.merge(left=res_df, right=mergeRes,
                          how='left', on=well_id_col)

        # Calculate thickness as % of total layer thickness (and clip to 1)
        res_df['TARG_THICK_PER'] = pd.DataFrame(np.round(res_df['TARG_THICK']/res_df['LAYER_THICK'], 3)).clip(upper=1)
        
        # Replace np.inf and np.nans with 0 (to ease interpolation)
        res_df['TARG_THICK_PER'] = res_df['TARG_THICK_PER'].where(res_df['TARG_THICK_PER'] != np.inf, other=0)
        res_df['TARG_THICK_PER'] = res_df['TARG_THICK_PER'].where(res_df['TARG_THICK_PER'] != np.nan, other=0)

        res_df["LAYER"] = layer  # Include the layer # as a separate gdf column
        res_df = res_df[[well_id_col, 'LATITUDE', 'LONGITUDE', 'LATITUDE_PROJ',
                         'LONGITUDE_PROJ', 'TOP', 'BOTTOM', 'TOP_ELEV',
                         'BOT_ELEV', 'SURFACE_ELEV', modelLyrTop_Col,
                         modelLyrBot_col, 'LAYER_THICK',
                         'TARG_THICK', 'TARG_THICK_PER', 'LAYER',
                         'geometry']].copy()
        
        # Convert back to GeoDataFrame
        res_df = gpd.GeoDataFrame(res_df,
                                  geometry='geometry',
                                  crs=gdf.crs)

        # Add layer outputs to list to be returned
        resdf_list.append(res_df)
        res_list.append(res)

        # Export data if specified
        if isinstance(export_dir, (pathlib.PurePath, str, os.PathLike)):
            export_dir = pathlib.Path(export_dir)
            if export_dir.is_dir():
                pass
            else:
                try:
                    os.mkdir(export_dir)
                except Exception:
                    print('Specified export directory does not exist and cannot be created. Function will continue to run, but data will not be exported.')

            # Format and build export filepath
            zFillDigs = len(str(len(layerList)))
            if str(outfile_prefix).endswith('_'):
                outfile_prefix = outfile_prefix[:-1]
            elif outfile_prefix is None:
                outfile_prefix = ''

            nowStr = str(datetime.datetime.today().date())+'_'+str(datetime.datetime.today().hour)+'-'+str(datetime.datetime.today().minute)+'-'+str(datetime.datetime.today().second)
            outPath = export_dir.joinpath(outfile_prefix+'_Lyr'+str(layer).zfill(zFillDigs)+'_'+nowStr+'.csv')
            outPathStr = outPath.as_posix()
            
            if return_all:
                res.to_csv(outPathStr, sep=',',
                           na_rep=np.nan, index_label='ID')
                res_df.to_csv(outPathStr, sep=',',
                              na_rep=np.nan, index_label='ID')
            else:
                res_df.to_csv(outPathStr, sep=',',
                              na_rep=np.nan, index_label='ID')

    if return_all:
        return res_list, resdf_list
    else:
        return resdf_list


# Interpolate layers to model model_grid
def layer_interp(points, model_grid, layers=None, interp_kind='nearest',
                 surface_grid=None, bedrock_grid=None,
                 layer_thick_grid=None, drift_thick_grid=None,
                 return_type='dataset', export_dir=None,
                 target_col='TARG_THICK_PER', layer_col='LAYER',
                 xcol=None, ycol=None, xcoord='x', ycoord='y',
                 log=False, verbose=False, **kwargs):

    """Function to interpolate results (by default TARG_THICK_PER,
    or Target thickness percent per layer).
    Converts dataframe points to gridded data.
    Results are saved to Model_Layer variable
    in output xarray.Dataset, by default
    
    This function uses the scipy.interpolate module for interpolation.

     Different interpolation methods may be used by specifying `interp_kind=`:
     * `'Nearest'`: Nearest neighbor (fastest).
     Uses scipy.interpolate.NearestNDInterpolator()
     * `'Linear'`: Linear interpolation
     Uses scipy.interpolate.LinearNDInterpolator()
     * `'Inter2d'`: Spline interpolation
     Uses scipy.interpolate.bisplrep()
     * `'CloughTocher`': Cubic interpolation using clough-tocher method
     Uses scipy.interpolate.CloughTocher2DInterpolator()
     * `'Radial basis function'`: Radial basis function
     Uses scipy.interpolate.RBFInterpolator()

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

    # Possible inputs (casefolded) for different interp methods
    nnList = ['nearest', 'nearest neighbor', 'nearestneighbor',
              'neighbor', 'nn', 'n']
    splineList = ['interp2d', 'interp2', 'interp', 'spline', 'spl', 'sp', 's']
    linList = ['linear', 'lin', 'l']
    ctList = ['clough tocher', 'clough', 'cloughtocher', 'ct', 'c']
    rbfList = ['rbf', 'radial basis', 'radial basis function', 'r', 'radial']

    # Potential future additions:
    #   k-nearest neighbors from scikit-learn?
    #   kriging? (from pykrige or maybe also from scikit-learn)

    # Check and format input model grid
    if not isinstance(model_grid, (xr.DataArray, xr.Dataset)):
        if pathlib.Path(model_grid).exists():
            model_grid = xr.open_dataarray(model_grid)

    # Extract model grid coordinates
    X = np.round(model_grid[xcoord].values, 10)
    Y = np.round(model_grid[ycoord].values, 10)

    # Get number of layers
    # If layers is not specified, use the length of the points list
    if layers is None and (type(points) is list or type(points) is dict):
        layers = len(points)

    if len(points) != layers:
        print('You have specified a different number of layers than what is iterable in the points argument. This may not work properly.')

    if verbose:
        print('\tInterpolating target lithology at each layer:')

    # Loop through each layer and interpolate (2d interpolation within layer)
    daDict = {}
    for lyr in range(1, layers+1):

        # Get points for each layer
        if type(points) is list or type(points) is dict:
            pts = points[lyr-1]
        else:
            pts = points

        # Get x coordinates for each well point
        if xcol is None:
            if 'geometry' in pts.columns:
                dataX = pts['geometry'].x
            else:
                print('xcol not specified and geometry column not detected (points is not/does not contain geopandas.GeoDataFrame)')
                return
        else:
            dataX = pts[xcol]
        
        # Get x coordinates for each well point
        if ycol is None:
            if 'geometry' in pts.columns:
                dataY = pts['geometry'].y
            else:
                print('ycol not specified and geometry column not detected (points is not/does not contain geopandas.GeoDataFrame)')
                return
        else:
            dataY = pts[ycol]

        # Get 'z' coordinates (interpolated value) for each well point
        interpVal = pts[target_col]

        # Drop points without x, y, or z data (and maintain consistency)
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

        # Carry out interpolation
        # Nearest neighbor
        if interp_kind.lower() in nnList:
            interpType = 'Nearest Neighbor'
            X, Y = np.meshgrid(X, Y, sparse=True)  # 2D Grid for interpolation
            dataPoints = np.array(list(zip(dataX, dataY)))
            interp = interpolate.NearestNDInterpolator(dataPoints, interpVal,
                                                       **kwargs)
            Z = interp(X, Y)

        # Linear
        elif interp_kind.lower() in linList:
            interpType = 'Linear'
            dataPoints = np.array(list(zip(dataX, dataY)))
            interp = interpolate.LinearNDInterpolator(dataPoints, interpVal,
                                                      **kwargs)
            X, Y = np.meshgrid(X, Y, sparse=True)  # 2D Grid for interpolation
            Z = interp(X, Y)

        # Clough-toucher (cubic)
        elif interp_kind.lower() in ctList:
            interpType = 'Clough-Toucher'
            X, Y = np.meshgrid(X, Y, sparse=True)  # 2D Grid for interpolation
            if 'tol' not in kwargs:
                kwargs['tol'] = 1e10
            interp = interpolate.CloughTocher2DInterpolator(list(zip(dataX, dataY)), interpVal, **kwargs)
            Z = interp(X, Y)

        # Radial basis function
        elif interp_kind.lower() in rbfList:
            interpType = 'Radial Basis'
            dataXY = np.column_stack((dataX, dataY))
            interp = interpolate.RBFInterpolator(dataXY, interpVal, **kwargs)
            print("Radial Basis Function does not work well with many well-based datasets. Consider instead specifying 'nearest', 'linear', 'spline', or 'clough tocher' for interpolation interp_kind.")
            Z = interp(np.column_stack((X.ravel(), Y.ravel()))).reshape(X.shape)

        # Spline
        elif interp_kind.lower() in splineList:
            interpType = 'Spline Interpolation'
            Z = interpolate.bisplrep(dataX, dataY, interpVal, **kwargs)

        # Nearest neighbor by default if otherwise not specified
        else:
            if verbose:
                print(f'Specified interpolation (interp_kind={interp_kind}) not recognized, using nearest neighbor.')
                print("\tinterp_kind should be one of: 'nearest', 'linear', 'interp2d', 'cloughtocher', or 'spline'")

            interpType = 'Nearest Neighbor'
            X, Y = np.meshgrid(X, Y, sparse=True)  # 2D Grid for interpolation
            interp = interpolate.NearestNDInterpolator(list(zip(dataX, dataY)),
                                                       interpVal, **kwargs)
            Z = interp(X, Y)

        # Create new datarray with new data values, else same as model_grid
        interp_grid = xr.DataArray(
                    data=Z,
                    dims=model_grid.dims,
                    coords=model_grid.coords)

        # Drop if only a single coordinate in "band" dimension
        if 'band' in interp_grid.coords:
            interp_grid = interp_grid.drop_vars('band')

        # Clip to 0-1 if percentage
        if target_col == 'TARG_THICK_PER':
            interp_grid = interp_grid.clip(min=0, max=1, keep_attrs=True)

        # Add coordinate in layer dimension for current layer
        interp_grid = interp_grid.expand_dims(dim='Layer')
        interp_grid = interp_grid.assign_coords(Layer=[lyr])

        # Delete to reduce memory usage (?)
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

    # Determine whether to export xarray DataArray or Dataset
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
        w4h.export_grids(grid_data=interp_data, out_path=export_dir, file_id='', filetype='tif', variable_sep=True, date_stamp=True)
        print('Exported to {}'.format(export_dir))

    return interp_data


def natural_neighbor_interp(points, model_grid=None,
                            well_id='API_NUMBER', interp_value='ELEVATION',
                            xcoord='LONGITUDE', ycoord='LATITUDE', elev="ELEVATION",
                            show_points=False,
                            show_plot=False, show_adjacent_regions=False, show_voronoi=False,
                            time_segments=False, layers=None, ):
    time0 = datetime.datetime.now()

    uniqueWellDF = points.drop_duplicates(subset=well_id)[[well_id, xcoord, ycoord, elev]].reset_index(drop=True)

    if show_points:
        uniqueWellDF.plot(xcoord, ycoord, kind='scatter')

    wellPtList = []
    wellValueList = []
    for i, well in uniqueWellDF.iterrows():
        wellPtList.append([well[xcoord], well[ycoord]])
        wellValueList.append(float(well[interp_value]))

    vClass = Voronoi(points=points)
    if time_segments:
        thisTime = lastTime = datetime.datetime.now()
        print(f"Voronoi time:  {thisTime-time0}")

    if show_voronoi:
        voronoi_plot_2d(vClass)

    # Get vertices as variable
    verts = vClass.vertices

    # Iterate through all coordinates for interpolation
    # CODE HERE FOR ITERATING THROUGH COORDS?
    newPoint = Point([-90.1, 38.7])

    def _vectorized_nat_neighbor(x, y):
        newPoint = Point([x, y])

        # Find region containing the "newPoint" (point for interpolation)
        polygonList = []
        polyRegInd = []
        containingRegion = False
        for i, region in enumerate(vClass.regions):
            # Create (valid) shapely polygon for each region using vertices
            vertList = [verts[v].tolist() for v in region if v != -1]

            if len(vertList) > 0 and vertList[0] != vertList[-1]:
                vertList.append(vertList[0])

            if len(vertList) == 0:
                continue
            currPoly = Polygon(vertList)

            # If polygon contains the point, record and break loop
            if currPoly.is_valid and currPoly.contains(newPoint) and containingRegion is False:
                containingRegion = currPoly
                containingRegionIndex = i
            polygonList.append(currPoly)
            polyRegInd.append(i)

        if time_segments:
            thisTime = datetime.datetime.now()
            print(f"Make Initial Polygons time:  {thisTime-lastTime}")
            lastTime = datetime.datetime.now()

        # Get all regions that touch containing region
        # UPDATE THIS WITH GEOPANDAS FOR EFFICIENCY!!!!!!
        regionAndAdjacents = [containingRegion]
        regionIndices = [containingRegionIndex]
        for i, region in enumerate(polygonList):
            touchCondition = containingRegion.touches(region)

            if touchCondition:
                regionAndAdjacents.append(region)
                regionIndices.append(polyRegInd[i])

        # Now, get second layer of regions
        # UPDATE THIS WITH GEOPANDAS FOR EFFICIENCY!!!!!!
        regionAndAdjacents2 = regionAndAdjacents.copy()
        for i, region in enumerate(regionAndAdjacents2):
            for j, regInd in enumerate(vClass.point_region):
                regionVerts = [verts[v].tolist() for v in vClass.regions[regInd]]
                currPoly = Polygon(regionVerts)
                if region.touches(currPoly) and regInd not in regionIndices:
                    regionAndAdjacents.append(currPoly)
                    regionIndices.append(int(regInd))

        if time_segments:
            thisTime = datetime.datetime.now()
            print(f"Get all adjacent and semi-adjacent Polygons:  {thisTime-lastTime}")
            lastTime = datetime.datetime.now()

        # Get points and values for interpolation
        nearPoints = [newPoint.coords[0]]
        nearValues = [None]
        vPR = list(vClass.point_region)
        for ind in regionIndices:
            nearValues.append(wellValueList[vPR.index(ind)])
            nearPoints.append(vClass.points[vPR.index(ind)].tolist())
        regGS = gpd.GeoDataFrame(zip(regionIndices, regionAndAdjacents), columns=['RegionIndices', 'geometry'], crs=4326)

        if time_segments:
            thisTime = datetime.datetime.now()
            print(f"Extract points for interpolation:  {thisTime-lastTime}")
            lastTime = datetime.datetime.now()

        # Plot regions and points
        if show_adjacent_regions:
            x, y = zip(*nearPoints)
            fig, ax = plt.subplots(figsize=(10,10))
            ax.scatter(x, y, c=nearValues)

            regGS.plot(ax=ax, facecolor='#00000000', edgecolor='k')
            ax.scatter(x=newPoint.xy[0], y=newPoint.xy[1],c='k', marker='+')
            [ax.text(x=x[i], y=y[i], s=f"{nv:.0f}") for i, nv in enumerate(nearValues) if i!=0]
            plt.show()
            if time_segments:
                thisTime = datetime.datetime.now()
                print(f"Plot adjacent regions:  {thisTime-lastTime}")
                lastTime = datetime.datetime.now()

        # Create new voronoi only with points with regions nearby to new point (to reduce computational cost)
        nearbyVor = Voronoi(points=nearPoints)
        if time_segments:
            thisTime = datetime.datetime.now()
            print(f"Generate near voronoi:  {thisTime-lastTime}")
            lastTime = datetime.datetime.now()

        # Go through just the nearby regions, check validitiy, and get info
        # (Might be able to just access with indices and save time?)
        polygons = []
        for i, region in enumerate(nearbyVor.regions):
            # Create (valid) shapely polygon for each region using vertices
            vertList = [verts[v].tolist() for v in region if v!=-1]

            if len(vertList)>0 and vertList[0] != vertList[-1]:
                vertList.append(vertList[0])

            if len(vertList)==0:
                continue
            currPoly = Polygon(vertList)

            # If polygon contains the point, record and break loop
            if currPoly.is_valid and currPoly.contains(newPoint) and containingRegion is False:
                containingRegion = currPoly
                containingRegionIndex = i
            polygons.append(currPoly)
            polyRegInd.append(i)

        if time_segments:
            thisTime = datetime.datetime.now()
            print(f"Get nearby polygons:  {thisTime-lastTime}")
            lastTime = datetime.datetime.now()

        # Get the "New" polygon (the one created for the interpolation point of interest)
        polyOfInterest = Polygon([nearbyVor.vertices[v] for v in nearbyVor.regions[nearbyVor.point_region[0]]])

        # Plot regions nearby point of interest, and overlap of new Voronoi and old
        if show_plot:
            fig3, ax3 = plt.subplots()
            regGS.plot(facecolor="#00000000", edgecolor='k', ax=ax3)
            x, y = zip(*nearPoints)
            ax3.scatter(x, y)
            gpd.GeoSeries([polyOfInterest], crs=4326).plot(facecolor="#00000000", 
                                                        edgecolor='r', linewidths=5, ax=ax3)
            if time_segments:
                thisTime = datetime.datetime.now()
                print(f"Plot nearby polys and overlap:  {thisTime-lastTime}")
                lastTime = datetime.datetime.now()

        # Check which regions overlap and their areas
        intersectionAreas = []
        intersectionValues = []
        for row, reg in regGS.iterrows():
            if reg['geometry'].intersects(polyOfInterest):
                intersectionAreas.append(polyOfInterest.intersection(reg['geometry']).area)
                intersectionValues.append(float(wellValueList[vPR.index(reg['RegionIndices'])]))            

        # Calculate weights
        weights = np.array(intersectionAreas)/np.sum(intersectionAreas)

        if time_segments:
            thisTime = datetime.datetime.now()
            print(f"Calculate weights:  {thisTime-lastTime}")
            lastTime = datetime.datetime.now()

        # Calculate and return natural neighbor interpolated value
        return float(np.sum(weights * intersectionValues))

    def coord_function_vec(lat, lon):
        return np.sin(np.radians(lat)) + np.cos(np.radians(lon))

    # Apply function vectorized over coordinates
    result_vec = xr.apply_ufunc(
        _vectorized_nat_neighbor,
        model_grid["x"].values[:, np.newaxis],  # broadcast along lon
        model_grid["y"].values[np.newaxis, :],  # broadcast along lat
        vectorize=True
        )

    return xr.DataArray(result_vec,
                        coords={"x": model_grid.y, "y": model_grid.x},
                        dims=("y", "x"))


# Optional, combine dataset
def combine_dataset(layer_dataset, surface_elev, bedrock_elev, layer_thick, log=False):
    """Function to combine xarray datasets or datarrays into
    a single xarray.Dataset.
    Useful to add surface, bedrock, layer thick, and layer datasets
    all into one variable, for pickling or netcdf export, for example.

    Parameters
    ----------
    layer_dataset : xr.DataArray
        DataArray contining all the interpolated layer information.
    surface_elev : xr.DataArray
        DataArray containing surface elevation data
    bedrock_elev : xr.DataArray
        DataArray containing bedrock elevation data
    layer_thick : xr.DataArray
        DataArray containing layer thickness at each point in the model grid
    log : bool, default = False
        Whether to log inputs and outputs to log file.

    Returns
    -------
    xarray.Dataset or xarray.DataArray
        Dataset with all arrays set to different variables within dataset.
        Or simply DataArray, if specified.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    daDict = {}
    daDict['Layers'] = layer_dataset
    daDict['Surface_Elev'] = surface_elev
    daDict['Bedrock_Elev'] = bedrock_elev
    daDict['Layer_Thickness'] = layer_thick

    combined_dataset = xr.Dataset(daDict)

    return combined_dataset
