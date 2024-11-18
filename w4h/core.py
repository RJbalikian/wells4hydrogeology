"""The Core module contains core functions of the package used in other modules or as primary functions in the package. 
This includes the main run() function that allows rapid data analysis, a function to retrieve sample data,
and functions that are used throughout the package for logging and printing verbose outputs."""

import datetime
import inspect
import json
import logging
import pathlib
import pkg_resources
import zipfile

import geopandas as gpd
import numpy as np
import pandas as pd
import rioxarray as rxr
from shapely import wkt
import xarray as xarray

from shapely.geometry import Point

import w4h

# Main function to run model all at once
def run(well_data,
        surf_elev_grid,
        bedrock_elev_grid,
        model_grid=None,
        metadata=None,
        keep_all_cols=True,
        layers = 9,
        description_col='FORMATION', top_col='TOP', bottom_col='BOTTOM', depth_type='depth',
        study_area=None, xcol='LONGITUDE', ycol='LATITUDE', zcol='SURFACE_ELEV', well_id_col='API_NUMBER',
        lith_dict=None, lith_dict_start=None, lith_dict_wildcard=None,
        target_dict=None,
        target_name='',
        include_elevation_grids=True,
        include_elevation_coordinates=True,
        export_dir=None,
        verbose=False,
        log=False,
        **kw_params):
    """Function to run entire process with one line of code. 
    
    NOTE: verbose and log are boolean parameters used for most of the functions. verbose=True prints information to terminal, log=True logs to a file in the log_dir, which defaults to the export_dir

    Parameters
    ----------
    well_data : str or pathlib.Path obj
        Filepath to file or directory containing well data.
    surf_elev_grid : str or pathlib.Path object
        _description_
    bedrock_elev_grid : str or pathlib.Path object
        _description_
    model_grid : str or pathlib.Path object, or model grid parameters (see model_grid function)
        _description_        
    metadata : str or pathlib.Path object, or None, default=None
        Filepath to file or directory containing well metadata, such as location and elevation. If None, will check if well_data is a directory, and if so, will use metadata_filename to search in same directory.
    keep_all_cols : bool, default=True
        Whether to keep all columns of the input dataframes/files. If True, no columns are excluded. If False, only keeps necessary columns.
    layers : int, default = 9
        The number of layers in the model grid
    description_col : str, default = 'FORMATION'
        Name of column containing geologic descriptions of the well interval. This column should be in well_data.
    top_col : str, default = 'TOP'
        Name of column containing depth/elevation at top of well interval. This column should be in well_data.
    bottom_col : str, default = 'BOTTOM'
        Name of column containing depth/elevation at bottom of well interval. This column should be in well_data.    
    depth_type : str, default = 'depth'
        Whether values top_col or bottom_col refer to depth or elevation.
    study_area : str or pathlib.Path object, or geopandas.GeoDataFrame
        _description_
    xcol : str, default = 'LONGITUDE' 
        Name of column containing x coordinates. This column should be in metadata unless metadata is not read, then it should be in well_data.
    ycol : str, default = 'LATITUDE'
        Name of column containing y coordinates. This column should be in metadata unless metadata is not read, then it should be in well_data.
    zcol : str, default = 'ELEVATION' 
        Name of column containing z coordinates. This column should be in metadata unless metadata is not read, then it should be in well_data.
    lith_dict : str or pathlib.Path object, or pandas.DataFrame
        _description_
    lith_dict_start : str or pathlib.Path object, or pandas.DataFrame
        _description_
    lith_dict_wildcard : str or pathlib.Path object, or pandas.DataFrame
        _description_
    target_dict : str or pathlib.Path object, or pandas.DataFrame
        _description_
    target_name : str, default = 'CoarseFine'
        Name of target of interest, to be used on exported files
    include_elevation_grid : bool, default = True
        Whether to include the elevation grids (surface, bedrock, and derived total "drift" thickness and layer thickness)
    include_elevation_coordinates : bool, default = True
        Whether to include the elevation coordinates for each grid point at all layers in the output.
        If True, these are saved as unindexed coordinates.
    export_dir : str or pathlib.Path object, default = None
        Directory to export output files
    verbose : bool, default = False
        Whether to print updates/results
    log : bool, default = False
        Whether to send parameters and outputs to log file, to be saved in export_dir, or the same directory as well_data if export_dir not defined.
    **kw_params
        Keyword parameters used by any of the functions throughout the process. See list of functions above, and the API documentation for their possible parameters
    """

    if verbose:
        verbose_print(run, locals())

    #Get data (files or otherwise)
    file_setup_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.file_setup).parameters.keys()}
    
    #Check how well_data and metadata were defined
    if isinstance(well_data, pathlib.PurePath) or isinstance(well_data, str):
        #Convert well_data to pathlib.Path if not already
        if isinstance(well_data, str):
            well_data = pathlib.Path(well_data)

        if metadata is None:
            if well_data.is_dir():
                # If the two files are supposed to be in the same directory (or just want well_data found)
                well_dataPath, metadataPath = w4h.file_setup(well_data=well_data, verbose=verbose, log=log, **file_setup_kwargs)             
            elif well_data.exists():
                # If well_data is a file, and metadata is not used
                well_dataPath, _ = w4h.file_setup(well_data=well_data, verbose=verbose, log=log, **file_setup_kwargs)             
                metadataPath = None
            else:
                #Need for well_data to exist at the very least
                raise IOError('well_data file does not exist:{}'.format(well_data))
        elif isinstance(metadata, pathlib.PurePath) or isinstance(metadata, str):
            #Metdata has specifically been specified by a filepath
            if isinstance(metadata, str):
                metadata = pathlib.Path(metadata)
            well_dataPath, metadataPath = w4h.file_setup(well_data=well_data, metadata=metadata, **file_setup_kwargs)                
        else:
            if isinstance(metadata, (pd.DataFrame, gpd.GeoDataFrame)):
                well_dataPath, _ = w4h.file_setup(well_data=well_data, verbose=verbose, log=log, **file_setup_kwargs)             
                metadataPath = metadata
            elif metadata is None:
                well_dataPath, _ = w4h.file_setup(well_data=well_data, verbose=verbose, log=log, **file_setup_kwargs)             

    elif isinstance(well_data, (pd.DataFrame, gpd.GeoDataFrame)):
        if isinstance(metadata, (pd.DataFrame, gpd.GeoDataFrame)):
            well_dataPath = well_data
            metadataPath = metadata
        elif isinstance(metadata, pathlib.PurePath) or isinstance(metadata, str):
            _, metadataPath = w4h.file_setup(well_data=metadata, metadata=metadata, verbose=verbose, log=log, **file_setup_kwargs)
            well_dataPath = well_data
        else:
            print('ERROR: metadata must be a string filepath, a pathlib.Path object, or pandas.DataFrame')
    else:
        print('ERROR: well_data must be a string filepath, a pathlib.Path object, or pandas.DataFrame')

    if not export_dir:
        if export_dir is False or export_dir is None:
            if verbose:
                print("\tData will not be exported")
            pass
        else:
            nowTime = datetime.datetime.now()
            nowTime = str(nowTime).replace(':', '-').replace(' ','_').split('.')[0]
            nowTimeStr = '_'+str(nowTime)
            outDir = 'Output_'+nowTimeStr
            if isinstance(well_dataPath, pd.DataFrame) or isinstance(well_dataPath, gpd.GeoDataFrame):
                export_dir = pathlib.Path(outDir)
            elif isinstance(well_dataPath, pathlib.PurePath):
                if well_dataPath.is_dir():
                    export_dir = well_dataPath.joinpath(outDir)
                else:
                    export_dir = well_dataPath.parent.joinpath(outDir)
            else:
                raise IOError('export_dir should be explicitly defined if well_data is not a filepath')

            if not export_dir.exists():
                try:
                    export_dir.mkdir()
                except Exception:
                    print('Export Directory not created')

    #Get pandas dataframes from input
    read_raw_txt_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.read_raw_csv).parameters.keys()}
    #if 'data_dtypes' not in read_raw_txt_kwargs.keys():
    #    read_raw_txt_kwargs['data_dtypes'] = w4h.read_dict(w4h.get_resources()['well_data_dtypes'])
    #if 'metadata_dtypes' not in read_raw_txt_kwargs.keys():
    #    read_raw_txt_kwargs['metadata_dtypes'] = w4h.read_dict(w4h.get_resources()['metadata_dtypes'])
    well_data_IN, metadata_IN = w4h.read_raw_csv(data_filepath=well_dataPath, metadata_filepath=metadataPath, verbose=verbose, log=log, **read_raw_txt_kwargs)
    # Functions to read data into dataframes. Also excludes extraneous columns, and drops header data with no location information

    # Define data types (file will need to be udpated)
    well_data_DF = w4h.define_dtypes(undefined_df=well_data_IN, datatypes=w4h.get_resources()['well_data_dtypes'], verbose=verbose, log=log)
    metadata_DF = w4h.define_dtypes(undefined_df=metadata_IN, datatypes=w4h.get_resources()['metadata_dtypes'], verbose=verbose, log=log)

    if metadata_DF is None:
        well_data_xyz = well_data_DF
    else:
        merge_metadata_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.merge_metadata).parameters.keys()}
        well_data_xyz = w4h.merge_metadata(data_df=well_data_DF, header_df=metadata_DF, data_cols=None, header_cols=None, auto_pick_cols=False, drop_duplicate_cols=True, log=False, **merge_metadata_kwargs)

    #Convert well_data_xyz to have geometry
    coords2geometry_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.coords2geometry).parameters.keys()}
    well_data_xyz = w4h.coords2geometry(df_no_geometry=well_data_xyz, xcol=xcol, ycol=ycol, zcol=zcol, verbose=verbose, log=log, **coords2geometry_kwargs)

    #Get Study area
    read_study_area_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.read_study_area).parameters.keys()}
    if study_area is None:
        studyAreaIN = None
        use_study_area = False
    else:
        studyAreaIN = w4h.read_study_area(study_area=study_area, log=log, **read_study_area_kwargs)
        use_study_area = True

    clip_gdf2study_area_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.clip_gdf2study_area).parameters.keys()}
    well_data_xyz = w4h.clip_gdf2study_area(study_area=studyAreaIN, gdf=well_data_xyz,  verbose=verbose, log=log,**clip_gdf2study_area_kwargs)
    #Get surfaces and grid(s)
    read_grid_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.read_grid).parameters.keys()}

    modelGridPath = model_grid
    surfaceElevPath = surf_elev_grid
    bedrockElevPath = bedrock_elev_grid

    modelGrid = w4h.read_grid(grid_path=modelGridPath, grid_type='model', study_area=studyAreaIN, verbose=verbose, log=log, **read_grid_kwargs)
    surfaceElevGridIN = w4h.read_grid(grid_path=surfaceElevPath, grid_type='surface', study_area=studyAreaIN, verbose=verbose, log=log, **read_grid_kwargs)
    bedrockElevGridIN = w4h.read_grid(grid_path=bedrockElevPath, grid_type='bedrock', study_area=studyAreaIN, verbose=verbose, log=log, **read_grid_kwargs)

    #UPDATE: MAKE SURE CRS's all align ***
    #Add control points
    add_control_points_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.add_control_points).parameters.keys()}
    well_data_xyz = w4h.add_control_points(df_without_control=well_data_xyz, xcol=xcol, ycol=ycol, zcol=zcol, top_col=top_col, bottom_col=bottom_col, description_col=description_col, verbose=verbose, log=log, **add_control_points_kwargs)

    #Analyze Surface(s) and grid(s)
    bedrockGrid, surfaceGrid = w4h.align_rasters(grids_unaligned=[bedrockElevGridIN, surfaceElevGridIN], model_grid=modelGrid, no_data_val_grid=0, log=log)
    driftThickGrid, layerThickGrid = w4h.get_drift_thick(surface_elev=surfaceGrid, bedrock_elev=bedrockGrid, layers=layers, plot=verbose, log=log)
    
    well_data_xyz = w4h.sample_raster_points(raster=bedrockGrid, points_df=well_data_xyz, xcol=xcol, ycol=ycol, new_col='BEDROCK_ELEV', verbose=verbose, log=log)
    well_data_xyz = w4h.sample_raster_points(raster=surfaceGrid, points_df=well_data_xyz, xcol=xcol, ycol=ycol, new_col='SURFACE_ELEV', verbose=verbose, log=log)
    well_data_xyz['BEDROCK_DEPTH'] = well_data_xyz['SURFACE_ELEV'] - well_data_xyz['BEDROCK_ELEV']
    well_data_xyz['LAYER_THICK'] = well_data_xyz['BEDROCK_DEPTH'] / layers
    
    well_data_xyz = w4h.get_layer_depths(df_with_depths=well_data_xyz, layers=layers, log=log)

    #Clean up data
    well_data_xyz = w4h.remove_nonlocated(df_with_locations=well_data_xyz, log=log, verbose=verbose)
    well_data_xyz = w4h.remove_no_topo(df_with_topo=well_data_xyz, zcol=zcol, verbose=verbose, log=log)

    remove_no_depth_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.remove_no_depth).parameters.keys()}
    well_data_xyz = w4h.remove_no_depth(well_data_xyz, verbose=verbose, top_col=top_col, bottom_col=bottom_col, log=log, **remove_no_depth_kwargs) #Drop records with no depth information

    remove_bad_depth_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.remove_bad_depth).parameters.keys()}
    well_data_xyz = w4h.remove_bad_depth(well_data_xyz, verbose=verbose, top_col=top_col, bottom_col=bottom_col, depth_type=depth_type, log=log, **remove_bad_depth_kwargs)#Drop records with bad depth information (i.e., top depth > bottom depth) (Also calculates thickness of each record)

    remove_no_formation_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.remove_no_description).parameters.keys()}
    well_data_xyz = w4h.remove_no_description(well_data_xyz, description_col=description_col, verbose=verbose, log=log, **remove_no_formation_kwargs)

    #CLASSIFICATION
    #Read dictionary definitions and classify
    get_search_terms_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.get_search_terms).parameters.keys()}
    specTermsPATH, startTermsPATH, wildcardTermsPATH, = w4h.get_search_terms(spec_path=lith_dict, start_path=lith_dict_start, wildcard_path=lith_dict_wildcard, log=log, **get_search_terms_kwargs)
    read_dictionary_terms_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.read_dictionary_terms).parameters.keys()}
    if 'class_flag' in read_dictionary_terms_kwargs.keys():
        del read_dictionary_terms_kwargs['class_flag'] #This is specific to an invidiual dict terms file, so don't want to use for all
    specTerms = w4h.read_dictionary_terms(dict_file=specTermsPATH, log=log, **read_dictionary_terms_kwargs)
    startTerms = w4h.read_dictionary_terms(dict_file=startTermsPATH, log=log, **read_dictionary_terms_kwargs)
    wildcardTerms = w4h.read_dictionary_terms(dict_file=wildcardTermsPATH, log=log, **read_dictionary_terms_kwargs)

    # Clean up dictionary terms
    specTerms = specTerms.drop_duplicates(subset='DESCRIPTION')
    specTerms = specTerms.reset_index(drop=True)

    startTerms = startTerms.drop_duplicates(subset='DESCRIPTION')
    startTerms = startTerms.reset_index(drop=True)

    wildcardTerms = wildcardTerms.drop_duplicates(subset='DESCRIPTION')
    wildcardTerms = wildcardTerms.reset_index(drop=True)

    if verbose:
        noSpecTerms = specTerms.shape[0]
        noStartTerms = startTerms.shape[0]
        noWildcardTerms = wildcardTerms.shape[0]
        
        print('\tSearch terms to be used:')
        print(f'\t\t {noSpecTerms} exact match term/definition pairs')
        print(f'\t\t {noStartTerms} starting match term/definition pairs')
        print(f'\t\t {noWildcardTerms} wildcard match term/definition pairs')

    #CLASSIFICATIONS
    #Exact match classifications
    well_data_xyz = w4h.specific_define(well_data_xyz, terms_df=specTerms, description_col=description_col, verbose=verbose, log=log)
    
    #.startswith classifications
    if lith_dict_start is not None:
        classifedDF, searchDF = w4h.split_defined(well_data_xyz, verbose=verbose, log=log)
        searchDF = w4h.start_define(df=searchDF, terms_df=startTerms, description_col=description_col, verbose=verbose, log=log)
        well_data_xyz = w4h.remerge_data(classifieddf=classifedDF, searchdf=searchDF) #UPDATE: Needed? ***    

    #wildcard/any substring match classifications
    if lith_dict_wildcard is not None:
        classifedDF, searchDF = w4h.split_defined(well_data_xyz, verbose=verbose, log=log)
        searchDF = w4h.wildcard_define(df=searchDF, terms_df=wildcardTerms, description_col=description_col, verbose=verbose, log=log)
        well_data_xyz = w4h.remerge_data(classifieddf=classifedDF, searchdf=searchDF) #UPDATE: Needed? ***    

    #Depth classification
    classifedDF, searchDF = w4h.split_defined(well_data_xyz, verbose=verbose, log=log)
    searchDF = w4h.depth_define(df=searchDF, thresh=550, verbose=verbose, log=log)
    well_data_xyz = w4h.remerge_data(classifieddf=classifedDF, searchdf=searchDF) #UPDATE: Needed? ***
    
    #Fill unclassified data
    well_data_xyz = w4h.fill_unclassified(well_data_xyz, classification_col='CLASS_FLAG')

    #Add target interpratations
    read_lithologies_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.read_lithologies).parameters.keys()}
    targetInterpDF = w4h.read_lithologies(lith_file=target_dict, log=log, **read_lithologies_kwargs)
    well_data_xyz = w4h.merge_lithologies(well_data_df=well_data_xyz, targinterps_df=targetInterpDF, target_col='TARGET', target_class='bool')

    #Sort dataframe to prepare for next steps
    #well_data_xyz = w4h.sort_dataframe(df=well_data_xyz, sort_cols=['API_NUMBER','TOP'], remove_nans=True)
    well_data_xyz = well_data_xyz.sort_values(by=[well_id_col, top_col])
    well_data_xyz = well_data_xyz.reset_index(drop=True)
    
    # UPDATE: Option to remove nans?
    well_data_xyz = well_data_xyz[well_data_xyz["LITHOLOGY"].notnull()]

    layer_target_thick_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.layer_target_thick).parameters.keys()}
    if 'return_all' in layer_target_thick_kwargs.keys():
        del layer_target_thick_kwargs['return_all'] #This needs to be set to False, so we don't want it reading in twice

    resdf = w4h.layer_target_thick(df=well_data_xyz, layers=layers, return_all=False, export_dir=export_dir, depth_top_col=top_col, depth_bot_col=bottom_col, log=log, **layer_target_thick_kwargs)
    
    # bedrockGrid, surfaceGrid, driftThickGrid, layerThickGrid
    layer_interp_kwargs = {k: v for k, v in locals()['kw_params'].items() if k in inspect.signature(w4h.layer_interp).parameters.keys()}
    layers_data = w4h.layer_interp(points=resdf, model_grid=modelGrid, layers=layers, verbose=verbose, log=log, **layer_interp_kwargs)

    if include_elevation_grids:
        # Add surface, bedrock, and derived grids
        layers_data['Surface_Elevation'] = surfaceGrid
        layers_data['Bedrock_Elevation'] = bedrockGrid
        layers_data['Bedrock_Depth'] = driftThickGrid
        layers_data['Layer_Thickness'] = layerThickGrid

    if include_elevation_coordinates:
        # Add each layer's elevation as an unindexed coordinate
        layerElevs = []
        for i in range(1, layers+1):
            layerElevs.append((layers_data['Surface_Elevation'] - (layers_data['Layer_Thickness']*i)).values)
        layerElevs = np.array(layerElevs)
        layers_data = layers_data.assign_coords(layer_elevs=(['Layer', "y", "x"], layerElevs))

    # Calculate current time for export string
    nowTime = datetime.datetime.now()
    nowTime = str(nowTime).replace(':', '-').replace(' ', '_').split('.')[0]
    nowTimeStr = '_'+str(nowTime)

    #THIS MAY BE REPEAT OF LAST LINES OF layer_interp()
    w4h.export_grids(grid_data=layers_data, out_path=export_dir, file_id=target_name,filetype='tif', variable_sep=True, date_stamp=True, verbose=verbose, log=log)

    return resdf, layers_data


# Function to update docstring for run function, used in __init__ file 
def _run_docstring():
    nl = '\n\t'
    functionList = [w4h.file_setup, w4h.read_raw_csv, w4h.define_dtypes, w4h.merge_metadata, w4h.coords2geometry,
                    w4h.read_study_area, w4h.clip_gdf2study_area, w4h.read_grid, w4h.add_control_points,
                    w4h.remove_nonlocated, w4h.remove_no_topo, w4h.remove_no_depth, w4h.remove_bad_depth, w4h.remove_no_description,
                    w4h.get_search_terms, w4h.read_dictionary_terms, w4h.specific_define, 
                    w4h.split_defined, w4h.start_define, w4h.wildcard_define, w4h.remerge_data, w4h.fill_unclassified,
                    w4h.read_lithologies, w4h.merge_lithologies, 
                    w4h.align_rasters, w4h.get_drift_thick, w4h.sample_raster_points, w4h.get_layer_depths, w4h.layer_target_thick,
                    w4h.layer_interp, w4h.export_grids]

    funcStrList = []
    funcParams = []
    funcDefaults = []
    prevOutputList = ['df', 'filepath', 'study_area']
    requiredList = []
    for func in functionList:
        parameters = inspect.signature(func).parameters
        defaults = [param.default for param in list(zip(*parameters.items()))[1]]
        parameters = list(zip(*parameters.items()))[0]

        for i, d in enumerate(defaults):
            if 'kwargs' in parameters[i]:
                defaults[i] = {}
            elif d is inspect._empty:
                if func.__name__ == 'read_study_area' and parameters[i] == 'study_area':
                    defaults[i] = "None <but defaults to w4h.resources()['study_area']>"
                elif any(o in parameters[i] for o in prevOutputList):
                    defaults[i] = '<output of previous function>'
                else:
                    defaults[i] = '<no default>'

        firstLine = f"\n\n**{func.__name__}**"
        followingLines = ''
        for i, param in enumerate(parameters):
            followingLines += f"\n\n\t{param}".ljust(25)
            if isinstance(defaults[i], str) and defaults[i]!='<positional, no default>':
                followingLines += f"| default = '{defaults[i]}'"
            else:
                followingLines += f"| default = {defaults[i]}"

        #funcDefaults.append(['<positional, no default>' if d is inspect._empty else d for d in defaults])
        #funcParams.append(list(zip(*parameters.items()))[0])

        funcString = firstLine + followingLines
        funcStrList.append(funcString)

    run_docstring = f"""
    w4h.run() is a function that runs the intended workflow of the wells4hydrogeology (w4h) package.
    This means that it runs several constituent functions. The workflow that this follows is provided in the package wiki.
    It accepts the parameters of the constituent functions. To see a list of these functions and parameters, use `help(w4h.run)`.

    The following functions used in w4h.run() are listed below, along with their parameters and default values for those parameters. 
    See the documentation for the each of the individual functions for more information on a specific parameter:

    {nl.join(funcStrList)}"

    """
    return run_docstring


# Function for logging (experimental)
log_filename = None  #initialize so variable exists but is None
def logger_function(logtocommence, parameters, func_name):
    """Function to log other functions, to be called from within other functions

    Parameters
    ----------
    logtocommence : bool
        Whether to perform logging steps
    parameters : dict
        Dictionary containing parameters and their values, from function
    func_name : str
        Name of function within which this is called
    """
    if logtocommence:
        global log_filename
        #log parameter should be false by default on all. If true, will show up in kwargs
        
        #Get the log parameter value
        if 'log' in parameters.keys():
            log_file = parameters.pop('log', None)
        else:
            #If it wasn't set, default to None
            log_file = None
        
        #Get currenet time and setup format for log messages
        curr_time = datetime.datetime.now()
        FORMAT = '%(asctime)s  %(message)s'

        #Check if we are starting a new logfile (only does this during run of file_setup() or (currently non-existent) new_logfile() functions)
        if log_file == True and (func_name == 'file_setup' or func_name == 'new_logfile'):

            #Get the log_dir variable set as a file_setup() parameter, or default to None if not specified
            out_dir = parameters.pop('log_dir', None)
            if out_dir is None:
                #If output directory not specified, default to the input directory
                out_dir = parameters['well_data']
            
            #Get the timestamp for the filename (this won't change, so represents the start of logging)
            timestamp = curr_time.strftime('%Y-%m-%d_%H-%M-%S')
            log_filename = pathlib.Path(out_dir).joinpath(f"log_{timestamp}.txt")
            if 'verbose' in parameters.keys():
                print('Logging data to', log_filename)

            #Set up logging stream using logging module
            logging.basicConfig(filename=log_filename, level=logging.INFO, format=FORMAT, filemode='w')

            #Log 
            logging.info(f"{func_name} CALLED WITH PARAMETERS:\n\t {parameters}")
        elif log_file == True:
            #Run this for functions that aren't setting up logging file
            if log_filename:
                #Get the log stream and log this function's call with parameters
                logging.basicConfig(filename=log_filename, level=logging.INFO, format=FORMAT)
                logging.info(f"{func_name} CALLED WITH PARAMETERS: \n\t{parameters}")
            else:
                #If log file has not already been set up, set it up
                timestamp = curr_time.strftime('%Y-%m-%d_%H-%M-%S')
                log_filename = f"log_{timestamp}.txt"

                #Now, get the log stream and log this function's call with parameters
                logging.basicConfig(filename=log_filename, level=logging.INFO, format=FORMAT)
                logging.info(f"{func_name} CALLED WITH PARAMETERS: \n\t{parameters}")
        else:
            #Don't log if log=False
            pass
    return


# Reusable function for consistently-formatted verbose printing output
def verbose_print(func, local_variables, exclude_params=[]):
    print_list = ['\n']
    sTime = datetime.datetime.now()
    print_list.append(f"{func.__name__}")
    print_list.append(f"\tStarted at {sTime}.")
    print_list.append(f"\tParameters:")
    for k, v in local_variables.items():
        if k in inspect.signature(func).parameters:
            if 'kwargs' in k:
                print_list.append(f"\t\t{k}")
                for kk, vv in local_variables[k].items():
                    print_list.append(f"\t\t\t{kk}={vv}")
            elif k in exclude_params:
                print_list.append(f"\t\t{k}=<input object>")
            else:
                print_list.append(f"\t\t{k}={v}")

    for line in print_list:
        print(line)
    return print_list


# Get filepaths for package resources in dictionary format
resource_dir = pathlib.Path(pkg_resources.resource_filename(__name__, 'resources/resources_home.txt')).parent
def get_resources(resource_type='filepaths', scope='local', verbose=False):
    """Function to get filepaths for resources included with package

    Parameters
    ----------
    resource_type : str, {'filepaths', 'data'}
        If filepaths, will return dictionary with filepaths to sample data. If data, returns dictionary with data objects.
    scope : str, {'local', 'statewide'}
        If 'local', will read in sample data for a local (around county sized) project. If 'state', will read in sample data for a statewide project (Illinois)
    verbose : bool, optional
        Whether to print results to terminal, by default False

    Returns
    -------
    resources_dict : dict
        Dictionary containing key, value pairs with filepaths to resources that may be of interest.
    """
    resources_dict = {}
    sample_data_dir = resource_dir.joinpath('sample_data')

    #Get sample data
    #Get lithology dictionaries' filepaths
    sample_dictionary_dir = sample_data_dir.joinpath('DictionaryTerms')
    resources_dict['LithologyDict_Exact'] = w4h.get_most_recent(dir=sample_dictionary_dir, glob_pattern='*DICTIONARY_SearchTerms*', verbose=verbose)
    resources_dict['LithologyDict_Start'] = w4h.get_most_recent(dir=sample_dictionary_dir, glob_pattern='*SearchTerms-Start*', verbose=verbose)
    resources_dict['LithologyDict_Wildcard'] = w4h.get_most_recent(dir=sample_dictionary_dir, glob_pattern='*SearchTerms-Wildcard*', verbose=verbose)

    #Get Lithology Interpretation filepaths
    lith_interp_dir = sample_data_dir.joinpath('LithologyInterpretations')
    resources_dict['LithInterps_FineCoarse'] = w4h.get_most_recent(dir=lith_interp_dir, glob_pattern='*FineCoarse*', verbose=verbose)
    resources_dict['LithInterps_Clay'] = w4h.get_most_recent(dir=lith_interp_dir, glob_pattern='*Clay*', verbose=verbose)
    resources_dict['LithInterps_Silt'] = w4h.get_most_recent(dir=lith_interp_dir, glob_pattern='*Silt*', verbose=verbose)    
    resources_dict['LithInterps_Sand'] = w4h.get_most_recent(dir=lith_interp_dir, glob_pattern='*Sand*', verbose=verbose)    
    resources_dict['LithInterps_Gravel'] = w4h.get_most_recent(dir=lith_interp_dir, glob_pattern='*Gravel*', verbose=verbose)    

    #Get other resource filepaths
    resources_dict['well_data_dtypes'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='*downholeDataTypes*', verbose=verbose)
    resources_dict['metadata_dtypes'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='*headerDataTypes*', verbose=verbose)
    resources_dict['ISWS_CRS'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='isws_crs.json', verbose=verbose)
    resources_dict['xyz_dtypes'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='xyzDataTypes.json', verbose=verbose)

    resources_dict['model_grid'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='grid_625_raster.tif', verbose=verbose)

    statewideSampleDir = sample_data_dir.joinpath('statewide_sample_data')
    statewideList = ['statewide', 'state', 'regional', 'region', 's', 'r']
    if scope.lower() in statewideList:
        resources_dict['well_data'] = statewideSampleDir.joinpath("IL_Statewide_WellData_XYz_2023-07-20_cleaned.zip")

        resources_dict['surf_elev'] = w4h.get_most_recent(dir=statewideSampleDir, glob_pattern='*IL_Statewide_Surface_Elev_ft_625ft_Lambert_GridAlign*', verbose=verbose)
        resources_dict['bedrock_elev'] = w4h.get_most_recent(dir=statewideSampleDir, glob_pattern='*IL_Statewide_Bedrock_Elev_2023_ft_625ft_Lambert_GridAlign*', verbose=verbose)
        resources_dict['study_area'] = w4h.get_most_recent(dir=statewideSampleDir, glob_pattern='*IL_Statewide_boundary*', verbose=verbose)
    else:
        resources_dict['study_area'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='*sample_studyArea*', verbose=verbose)
        resources_dict['surf_elev'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='*sample_surface_bedrock_lidarresampled100ft*', verbose=verbose)
        resources_dict['bedrock_elev'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='*LocalSample_Bedrock_elev_EStLGrimleyPhillips*', verbose=verbose)

        resources_dict['well_data'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='sample_well_data*', verbose=verbose)

    # Get data objects if specified
    dataObjList = ['data', 'objects', 'do', 'data objects', 'dataobjects']
    if resource_type.lower() in dataObjList:
        resources_dict['LithologyDict_Exact'] = pd.read_csv(resources_dict['LithologyDict_Exact'], 
                                                            dtype={"ID":int, "DESCRIPTION":str, "LITHOLOGY":str,
                                                            "COLOR":str, "CONSISTENCY":str, "MOD1":str, "MOD2":str,
                                                            "INTERPRETED":str, "COMPLETED":str, "ORIGIN_INDIANA":str},
                                                            index_col='ID')
        resources_dict['LithologyDict_Start'] = pd.read_csv(resources_dict['LithologyDict_Start'])
        resources_dict['LithologyDict_Wildcard'] = pd.read_csv(resources_dict['LithologyDict_Wildcard'])

        resources_dict['LithInterps_FineCoarse'] = pd.read_csv(resources_dict['LithInterps_FineCoarse'])
        resources_dict['LithInterps_Clay'] = pd.read_csv(resources_dict['LithInterps_Clay'])
        resources_dict['LithInterps_Silt'] = pd.read_csv(resources_dict['LithInterps_Silt'])
        resources_dict['LithInterps_Sand'] = pd.read_csv(resources_dict['LithInterps_Sand'])
        resources_dict['LithInterps_Gravel'] = pd.read_csv(resources_dict['LithInterps_Gravel'])

        
        with open(resources_dict['well_data_dtypes'], 'r', encoding='utf-8') as f:
            resources_dict['well_data_dtypes'] = json.load(f)

        with open(resources_dict['metadata_dtypes'], 'r', encoding='utf-8') as f:
            resources_dict['metadata_dtypes'] = json.load(f)            

        with open(resources_dict['ISWS_CRS'], 'r', encoding='utf-8') as f:
            resources_dict['ISWS_CRS'] = json.load(f)
        
        with open(resources_dict['xyz_dtypes'], 'r', encoding='utf-8') as f:
            resources_dict['xyz_dtypes'] = json.load(f)


        if scope.lower() in statewideList:
            sacrs = resources_dict['ISWS_CRS']
            with zipfile.ZipFile(resources_dict['well_data'].as_posix(), 'r') as archive:
                for file_name in archive.namelist():
                    with archive.open(file_name) as file:
                        if 'HEADER' in file_name:
                            metaDF = pd.read_csv(file)
                        else:
                            resources_dict['well_data'] = pd.read_csv(file)
            geometry = [Point(xy) for xy in zip(resources_dict['well_data']['LONGITUDE'], resources_dict['well_data']['LATITUDE'])]
            resources_dict['well_data'] = gpd.GeoDataFrame(resources_dict['well_data'], geometry=geometry, crs='EPSG:5070')
            
        else:
            sacrs = 'EPSG:5070'
            df = pd.read_csv(resources_dict['well_data'])
            df['geometry'] = df['geometry'].apply(wkt.loads)
            resources_dict['well_data'] = gpd.GeoDataFrame(df, geometry='geometry')


        resources_dict['study_area'] = gpd.read_file(resources_dict['study_area'], geometry='geometry', crs=sacrs)

        resources_dict['model_grid'] = rxr.open_rasterio(resources_dict['model_grid'])
        resources_dict['surf_elev'] = rxr.open_rasterio(resources_dict['surf_elev'])
        #resources_dict['surf_elev'] = resources_dict['surf_elev'].sel(band=1)
        resources_dict['bedrock_elev'] = rxr.open_rasterio(resources_dict['bedrock_elev'])
        #resources_dict['bedrock_elev'] = resources_dict['bedrock_elev'].sel(band=1)

    return resources_dict


# Only used for development purposes, check that parameters are unique
def __check_parameter_names(verbose=True):
    #Check parameters are unique
    import inspect
    import w4h
    import pandas as pd
    function_list = [w4h.file_setup,
                 w4h.read_raw_csv,
                 w4h.define_dtypes,
                 w4h.read_study_area,
                 w4h.read_grid,
                 w4h.add_control_points,
                 w4h.coords2geometry,
                 w4h.clip_gdf2study_area,
                 w4h.remove_nonlocated,
                 w4h.remove_no_topo,
                 w4h.remove_no_depth,
                 w4h.remove_bad_depth,
                 w4h.remove_no_description,
                 w4h.get_search_terms,
                 w4h.read_dictionary_terms,
                 w4h.specific_define,
                 w4h.start_define,
                 w4h.wildcard_define,
                 w4h.depth_define,
                 w4h.fill_unclassified,
                 w4h.read_lithologies,
                 w4h.merge_lithologies,
                 w4h.align_rasters,
                 w4h.get_drift_thick,
                 w4h.sample_raster_points,
                 w4h.get_layer_depths,
                 w4h.layer_target_thick,
                 w4h.layer_interp,
                 w4h.export_grids]
    
    paramDF = pd.DataFrame()
    for f in function_list:
        currParamList = inspect.getfullargspec(f)[0]
        fList = []
        for p in currParamList:
            fList.append(f.__name__)
        currParamDF = pd.DataFrame({'Function':fList, 'Parameter':currParamList})
        paramDF = pd.concat([paramDF, currParamDF])

    uniqueDF = paramDF.drop_duplicates(subset='Parameter').copy()

    for up in uniqueDF['Parameter']:
        if up != 'verbose' and up!='log':
            matchDF = paramDF[paramDF['Parameter']==up].copy()
            if verbose:
                if matchDF.shape[0] > 1:
                    print(matchDF)
    
    return paramDF
