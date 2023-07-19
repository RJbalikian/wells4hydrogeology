import datetime
import logging
import os
import pathlib

import pandas as pd

import w4h

import pkg_resources

def run(well_data, well_data_cols=None, 
        metadata=None, well_metadata_cols=None, 
        layers = 9,
        description_col='FORMATION', top_col='TOP', bottom_col='BOTTOM', depth_type='depth',
        study_area=None, xcol='LONGITUDE', ycol='LATITUDE', zcol='ELEVATION', well_id_col='API_NUMBER', output_crs='EPSG:4269',
        surf_elev_file=None, bedrock_elev_file=None, model_grid=None,
        lith_dict=None, lith_dict_start=None, lith_dict_wildcard=None,
        target_dict=None,
        target_name='CoarseFine',
        export_dir=None,
        verbose=False,
        log=False,
        **keyword_parameters):
    
    """Function to run entire process with one line of code. 
    
    NOTE: verbose and log are boolean parameters used for most of the functions. verbose=True prints information to terminal, log=True logs to a file in the log_dir, which defaults to the export_dir

    Parameters
    ----------
    well_data : str or pathlib.Path obj
        Filepath to file or directory containing well data.
    well_data_cols : List or list-like
        Columns to 
    metadata : str or pathlib.Path object
        Filepath to file or directory containing well metadata, such as location and elevation.
    well_metadata_cols : List or list-like
        _description_
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
    output_crs : crs definition accepted by pyproj, default = 'EPSG:4269'
        CRS to output all of the data into
    surf_elev_file : str or pathlib.Path object
        _description_
    bedrock_elev_file : str or pathlib.Path object
        _description_
    model_grid : str or pathlib.Path object, or model grid parameters (see model_grid function)
        _description_
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
    export_dir : str or pathlib.Path object, default = None
        Directory to export output files
    verbose : bool, default = False
        Whether to print updates/results
    log : bool, default = False
        Whether to send parameters and outputs to log file, to be saved in export_dir, or the same directory as well_data if export_dir not defined.
    **keyword_parameters
        Keyword parameters used by any of the functions throughout the process. See list of functions above, and the API documentation for their possible parameters
    """

    #Get data (files or otherwise)
    file_setup_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.file_setup.__code__.co_varnames}
    
    #Check how well_data and metadata were defined
    if isinstance(well_data, pathlib.PurePath) or isinstance(well_data, str):
        #Convert well_data to pathlib.Path if not already
        if isinstance(well_data, str):
            well_data = pathlib.Path(well_data)

        if metadata is None:
            if well_data.is_dir():
                #If the two files are supposed to be in the same directory (or just want well_data found)
                downholeDataPATH, headerDataPATH = w4h.file_setup(well_data=well_data, verbose=verbose, log=log, **file_setup_kwargs)             
            elif well_data.exists():
                #If well_data is a file, and metadata is not used
                downholeDataPATH, _ = w4h.file_setup(well_data=well_data, verbose=verbose, log=log, **file_setup_kwargs)             
                headerDataPATH = None
            else:
                #Need for well_data to exist at the very least
                IOError('well_data file does not exist:{}'.format(well_data))
        elif isinstance(metadata, pathlib.PurePath) or isinstance(metadata, str):
            #Metdata has specifically been specified by a filepath
            if isinstance(metadata, str):
                metadata = pathlib.Path(metadata)    
            downholeDataPATH, headerDataPATH = w4h.file_setup(well_data=well_data, metadata=metadata, **file_setup_kwargs)                
        else:
            if isinstance(metadata, pd.DataFrame):
                downholeDataPATH, _ = w4h.file_setup(well_data=well_data, verbose=verbose, log=log, **file_setup_kwargs)             
                headerDataPATH = metadata
            elif metadata is None:
                downholeDataPATH, _ = w4h.file_setup(well_data=well_data, verbose=verbose, log=log, **file_setup_kwargs)             

    elif isinstance(well_data, pd.DataFrame):
        if isinstance(metadata, pd.DataFrame):
            downholeDataPATH = well_data
            headerDataPATH = metadata
        elif isinstance(metadata, pathlib.PurePath) or isinstance(metadata, str):
            _, headerDataPATH = w4h.file_setup(well_data=metadata, metadata=metadata, verbose=verbose, log=log, **file_setup_kwargs)                
            downholeDataPATH = well_data
        else:
            print('ERROR: metadata must be a string filepath, a pathlib.Path object, or pandas.DataFrame')
    else:
        print('ERROR: well_data must be a string filepath, a pathlib.Path object, or pandas.DataFrame')

    #Get pandas dataframes from input
    read_raw_txt_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.read_raw_csv .__code__.co_varnames}
    well_data_IN, metadata_IN = w4h.read_raw_csv(data_filepath=downholeDataPATH, metadata_filepath=headerDataPATH, verbose=verbose, log=log, **read_raw_txt_kwargs) 
    #Functions to read data into dataframes. Also excludes extraneous columns, and drops header data with no location information

    #Define data types (file will need to be udpated)
    well_data_DF = w4h.define_dtypes(undefined_df=well_data_IN, datatypes='./resources/downholeDataTypes.txt', verbose=verbose, log=log)
    metadata_DF = w4h.define_dtypes(undefined_df=metadata_IN, datatypes='./resources/headerDataTypes.txt', verbose=verbose, log=log)

    if metadata_DF is None:
        well_data_xyz = well_data_DF
    else:
        merge_tables_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.merge_tables.__code__.co_varnames}
        well_data_xyz = w4h.merge_tables(data_df=well_data_DF, header_df=metadata_DF, data_cols=None, header_cols=None, auto_pick_cols=False, drop_duplicate_cols=True, log=False, **merge_tables_kwargs)

    #Convert well_data_xyz to have geometry
    coords2geometry_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.coords2geometry.__code__.co_varnames}
    well_data_xyz = w4h.coords2geometry(df_no_geometry=well_data_xyz, xcol=xcol, ycol=ycol, zcol=zcol, log=log, **coords2geometry_kwargs)

    #Get Study area
    read_study_area_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.read_study_area.__code__.co_varnames}
    if study_area is None:
        studyAreaIN = None
        use_study_area = False
    else:
        studyAreaIN = w4h.read_study_area(study_area_path=study_area, log=log, **read_study_area_kwargs)
        use_study_area = True

    clip_gdf2study_area_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.clip_gdf2study_area.__code__.co_varnames}
    well_data_xyz = w4h.clip_gdf2study_area(study_area=studyAreaIN, gdf=well_data_xyz, log=log, **clip_gdf2study_area_kwargs)

    #Get surfaces and grid(s)
    read_grid_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.read_grid.__code__.co_varnames}

    modelGridPath = model_grid
    surfaceElevPath = surf_elev_file
    bedrockElevPath = bedrock_elev_file
    #UPDATE: allow other types of model grid read ***
    modelGrid = w4h.read_grid(grid_path=modelGridPath, grid_type='model', study_area=studyAreaIN, verbose=verbose, log=log, **read_grid_kwargs)
    surfaceElevGridIN = w4h.read_grid(grid_path=surfaceElevPath, grid_type='surface', study_area=studyAreaIN, verbose=verbose, log=log, **read_grid_kwargs)
    bedrockElevGridIN = w4h.read_grid(grid_path=bedrockElevPath, grid_type='bedrock', study_area=studyAreaIN, verbose=verbose, log=log, **read_grid_kwargs)

    #UPDATE: MAKE SURE CRS's all align ***
    #Add control points
    #UPDATE: Code here for adding in control points ***
    add_control_points_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.add_control_points.__code__.co_varnames}
    well_data_xyz = w4h.add_control_points(df_without_control=well_data_xyz, xcol=xcol, ycol=ycol, zcol=zcol, top_col=top_col, bottom_col=bottom_col, description_col=description_col, verbose=verbose, log=log, **add_control_points_kwargs)
    
    #Clean up data
    well_data_xyz = w4h.remove_nonlocated(df_with_locations=well_data_xyz, log=log, verbose=verbose)
    well_data_xyz = w4h.remove_no_topo(df_with_topo=well_data_xyz, zcol=zcol, verbose=verbose, log=log)

    remove_no_depth_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.remove_no_depth.__code__.co_varnames}
    well_data_xyz = w4h.remove_no_depth(well_data_xyz, verbose=verbose, top_col=top_col, bottom_col=bottom_col, log=log, **remove_no_depth_kwargs) #Drop records with no depth information

    remove_bad_depth_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.remove_bad_depth.__code__.co_varnames}
    well_data_xyz = w4h.remove_bad_depth(well_data_xyz, verbose=verbose, top_col=top_col, bottom_col=bottom_col, depth_type=depth_type, log=log, **remove_bad_depth_kwargs)#Drop records with bad depth information (i.e., top depth > bottom depth) (Also calculates thickness of each record)

    remove_no_formation_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.remove_no_description.__code__.co_varnames}
    well_data_xyz = w4h.remove_no_description(well_data_xyz, description_col=description_col, verbose=verbose, log=log, **remove_no_formation_kwargs)

    #CLASSIFICATION
    #Read dictionary definitions and classify
    get_search_terms_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.get_search_terms.__code__.co_varnames}
    specTermsPATH, startTermsPATH, wildcardTermsPATH, = w4h.get_search_terms(spec_path=lith_dict, start_path=lith_dict_start, wildcard_path=lith_dict_wildcard, log=log, **get_search_terms_kwargs)
    read_dictionary_terms_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.read_dictionary_terms.__code__.co_varnames}
    specTerms = w4h.read_dictionary_terms(dict_file=specTermsPATH, log=log, **read_dictionary_terms_kwargs)
    startTerms = w4h.read_dictionary_terms(dict_file=startTermsPATH, log=log, **read_dictionary_terms_kwargs)
    wildcardTerms = w4h.read_dictionary_terms(dict_file=wildcardTermsPATH, log=log, **read_dictionary_terms_kwargs)

    #Clean up dictionary terms
    specTerms.drop_duplicates(subset='DESCRIPTION', inplace=True)
    specTerms.reset_index(inplace=True, drop=True)

    startTerms.drop_duplicates(subset='DESCRIPTION', inplace=True)
    startTerms.reset_index(inplace=True, drop=True)

    wildcardTerms.drop_duplicates(subset='DESCRIPTION', inplace=True)
    wildcardTerms.reset_index(inplace=True, drop=True)

    if verbose:
        print('Search terms to be used:')
        print('\t {} exact match term/definition pairs')
        print('\t {} starting match term/definition pairs')
        print('\t {} wildcard match term/definition pairs')

    #CLASSIFICATIONS
    #Exact match classifications
    well_data_xyz = w4h.specific_define(well_data_xyz, specTerms, description_col=description_col, verbose=verbose, log=log)
    
    #.startswith classifications
    classifedDF, searchDF = w4h.split_defined(well_data_xyz, verbose=verbose, log=log)
    searchDF = w4h.start_define(df=searchDF, terms_df=startTerms, description_col=description_col, verbose=verbose, log=log)
    well_data_xyz = w4h.remerge_data(classifieddf=classifedDF, searchdf=searchDF) #UPDATE: Needed? ***    

    #wildcard/any substring match classifications    
    classifedDF, searchDF = w4h.split_defined(well_data_xyz, verbose=verbose, log=log)
    searchDF = w4h.wildcard_define(df=searchDF, terms_df=wildcardTerms, description_col=description_col, verbose=verbose, log=log)
    well_data_xyz = w4h.remerge_data(classifieddf=classifedDF, searchdf=searchDF) #UPDATE: Needed? ***    

    #Depth classification
    classifedDF, searchDF = w4h.split_defined(well_data_xyz, verbose=verbose, log=log)
    searchDF = w4h.depth_define(searchDF, thresh=550, verbose=verbose, log=log)
    well_data_xyz = w4h.remerge_data(classifieddf=classifedDF, searchdf=searchDF) #UPDATE: Needed? ***
    
    #Fill unclassified data
    well_data_xyz = w4h.fill_unclassified(well_data_xyz, classification_col='CLASS_FLAG')
    
    #Add target interpratations
    read_lithologies_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.read_lithologies.__code__.co_varnames}
    targetInterpDF = w4h.read_lithologies(lith_file=target_dict, log=log, **read_lithologies_kwargs)
    well_data_xyz = w4h.merge_lithologies(df=well_data_xyz, targinterps_df=targetInterpDF, target_col='TARGET', target_class='bool')
    
    #Sort dataframe to prepare for next steps
    #well_data_xyz = w4h.sort_dataframe(df=well_data_xyz, sort_cols=['API_NUMBER','TOP'], remove_nans=True)
    well_data_xyz = well_data_xyz.sort_values(by=[well_id_col, top_col])
    well_data_xyz.reset_index(inplace=True, drop=True)
    #UPDATE: Option to remove nans?
    well_data_xyz = well_data_xyz[pd.notna(well_data_xyz["LITHOLOGY"])]

    #Analyze Surface(s) and grid(s)
    bedrockGrid, surfaceGrid = w4h.align_rasters(grids_unaligned=[bedrockElevGridIN, surfaceElevGridIN], modelgrid=modelGrid, no_data_val_grid=0, log=log)
    driftThickGrid, layerThickGrid = w4h.get_drift_thick(surface=surfaceGrid, bedrock=bedrockGrid, layers=layers, plot=verbose, log=log)
    #UPDATE: LAYER NAMES SO DON"T INCLUDE FT
    well_data_xyz = w4h.sample_raster_points(raster=bedrockGrid, points_df=well_data_xyz, xcol=xcol, ycol=ycol, new_col='BEDROCK_ELEV_FT', verbose=verbose, log=log)
    well_data_xyz = w4h.sample_raster_points(raster=surfaceGrid, points_df=well_data_xyz, xcol=xcol, ycol=ycol, new_col='SURFACE_ELEV_FT', verbose=verbose, log=log)
    well_data_xyz = w4h.sample_raster_points(raster=driftThickGrid, points_df=well_data_xyz, xcol=xcol, ycol=ycol, new_col='BEDROCK_DEPTH_FT', verbose=verbose, log=log)
    well_data_xyz = w4h.sample_raster_points(raster=layerThickGrid, points_df=well_data_xyz, xcol=xcol, ycol=ycol, new_col='LAYER_THICK_FT', verbose=verbose, log=log)
    well_data_xyz = w4h.get_layer_depths(df_with_depths=well_data_xyz, no_layers=layers, log=log)

    resdf = w4h.layer_target_thick(well_data_xyz, layers=9, return_all=False, outfile_prefix='CoarseFine', export_dir=export_dir, depth_top_col=top_col, depth_bot_col=bottom_col, log=log)
    
    layer_interp_kwargs = {k: v for k, v in locals()['keyword_parameters'].items() if k in w4h.layer_interp.__code__.co_varnames}
    layers_data = w4h.layer_interp(points=resdf, layers=9, grid=modelGrid, verbose=verbose, log=log, **layer_interp_kwargs)

    if export_dir is None:
        if well_data.is_dir():
            export_dir = well_data.joinpath('Output')
        else:
            export_dir = well_data.parent.joinpath('Outpath')
        
        if not export_dir.exists():
            try:
                export_dir.mkdir()
            except:
                pass

    w4h.export_grids(grid_data=layers_data, out_path=export_dir, file_id='',filetype='tif', variable_sep=True, date_stamp=True, verbose=verbose, log=log)
    #UPDATE: export points?
    return resdf, layers_data


log_filename=None #Set up so exists but is None
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

#Get filepaths for package resources in dictionary format
resource_dir = pathlib.Path(pkg_resources.resource_filename(__name__, 'resources/resources_home.txt')).parent
def get_resources(verbose=False):
    """Function to get filepaths for resources included with package

    Parameters
    ----------
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
    resources_dict['well_data_dtypes'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='downholeDataTypes.txt', verbose=verbose)
    resources_dict['metadata_dtypes'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='headerDataTypes.txt', verbose=verbose)
    resources_dict['ISWS_CRS'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='isws_crs.txt', verbose=verbose)
    resources_dict['xyz_dtypes'] = w4h.get_most_recent(dir=sample_data_dir, glob_pattern='xyzDataTypes.txt', verbose=verbose)

    return resources_dict

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
