import datetime
import logging
import os
import pathlib

import pandas as pd

import w4h

#log_filename = None

#Log data to file
"""def logger(func):
    def wrapper(*args, **kwargs):
        global log_filename
        #log parameter should be false by default on all. If true, will show up in kwargs
            #Is there a way to do this so all can be set at once?
        if 'log' in kwargs.keys():
            log_file = kwargs.pop('log', None)
        else:
            log_file = None
        if log_file == True and (func.__name__ == 'file_setup' or func.__name__ == 'new_logfile'):
            out_dir = kwargs.pop('out_dir', None)
            if out_dir is None:
                out_dir = kwargs['db_dir']
            timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
            log_filename = f"log_{timestamp}.txt"
            logging.basicConfig(filename=log_filename, level=logging.INFO)
        elif log_file == True:
            if log_filename:
                logging.basicConfig(filename=log_filename, level=logging.INFO)
            else:
                timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
                log_filename = f"log_{timestamp}.txt"
                logging.basicConfig(filename=log_filename, level=logging.INFO)
        else:
            pass
        result = func(*args, **kwargs)
        print('logged', func.__name__)
        print('fname', log_filename)
        return result
    return wrapper"""

log_filename=None

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
            #Is there a way to do this so all can be set at once?
        if 'log' in parameters.keys():
            log_file = parameters.pop('log', None)
        else:
            log_file = None
        
        curr_time = datetime.datetime.now()
        FORMAT = '%(asctime)s  %(message)s'
        if log_file == True and (func_name == 'file_setup' or func_name == 'new_logfile'):
            out_dir = parameters.pop('log_dir', None)
            if out_dir is None:
                out_dir = parameters['db_dir']
            timestamp = curr_time.strftime('%Y-%m-%d_%H-%M-%S')
            log_filename = pathlib.Path(out_dir).joinpath(f"log_{timestamp}.txt")
            print('Logging data to', log_filename)
            logging.basicConfig(filename=log_filename, level=logging.INFO, format=FORMAT, filemode='w')
            logging.info(f"Called {func_name} with args: {parameters}")
        elif log_file == True:
            if log_filename:
                logging.basicConfig(filename=log_filename, level=logging.INFO, format=FORMAT)
                logging.info(f"Called {func_name} with args: {parameters}")
            else:
                timestamp = curr_time.strftime('%Y-%m-%d_%H-%M-%S')
                log_filename = f"log_{timestamp}.txt"
                logging.basicConfig(filename=log_filename, level=logging.INFO, format=FORMAT)
                logging.info(f"Called {func_name} with args: {parameters}")
        else:
            pass
    return

def run(well_data, well_data_cols=None, 
        well_metadata=None, well_metadata_cols=None, 
        description_col='FORMATIoN', top_col='TOP', bottom_col='BOTTOM', depth_type='depth',
        xcol='LONGITUDE', ycol='LATITUDE', zcol='ELEVATION', idcol='API_NUMBER', output_crs='EPSG:4269',
        surf_elev_file=None, bedrock_elev_file=None, model_grid=None,
        lith_dict=None, lith_dict_start=None, lith_dict_wildcard=None,
        target_dict=None,
        study_area=None,
        export_dir=None,
        verbose=False,
        log=False,
        **keyword_parameters):
    
    """Function to run entire process with one line of code

    Parameters
    ----------
    well_data : str or pathlib.Path obj
        Filepath to file or directory containing well data.
    well_data_cols : List or list-like
        Columns to 
    well_metadata : str or pathlib.Path object
        Filepath to file or directory containing well metadata, such as location and elevation.
    well_metadata_cols : List or list-like
        _description_
    description_col : str, default = 'FORMATION'
        Name of column containing geologic descriptions of the well interval. This column should be in well_data.
    top_col : str, default = 'TOP'
        Name of column containing depth/elevation at top of well interval. This column should be in well_data.
    bottom_col : str, default = 'BOTTOM'
        Name of column containing depth/elevation at bottom of well interval. This column should be in well_data.    
    depth_type : str, default = 'depth'
        Whether values top_col or bottom_col refer to depth or elevation.
    xcol : str, default = 'LONGITUDE' 
        Name of column containing x coordinates. This column should be in well_metadata unless well_metadata is not read, then it should be in well_data.
    ycol : str, default = 'LATITUDE'
        Name of column containing y coordinates. This column should be in well_metadata unless well_metadata is not read, then it should be in well_data.
    zcol : str, default = 'ELEVATION' 
        Name of column containing z coordinates. This column should be in well_metadata unless well_metadata is not read, then it should be in well_data.
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
    study_area : str or pathlib.Path object, or geopandas.GeoDataFrame
        _description_
    """

    #Get important information
    todayDate, dateSuffix = w4h.get_current_date() 
    repoDir = pathlib.Path(os.getcwd()) #this will need to be updated for pypi packaging

    #Get data (files or otherwise)
    file_setup_kwargs = {k: v for k, v in locals().items() if k in w4h.file_setup.__code__.co_varnames}
    #Check how well_data was defined
    if isinstance(well_data, pathlib.PurePath) or isinstance(well_data, str):
        #Convert well_data to pathlib.Path if not already
        if isinstance(well_data, str):
            well_data = pathlib.Path(well_data)

        if isinstance(well_metadata, pathlib.PurePath) or isinstance(well_metadata, str):
            if isinstance(well_metadata, str):
                well_metadata = pathlib.Path(well_metadata)    
            downholeDataPATH, headerDataPATH = w4h.file_setup(well_data_path=well_data, metadata_path=well_metadata, verbose=verbose, log=log, **file_setup_kwargs)                
        else:
            if isinstance(well_metadata, pd.DataFrame):
                downholeDataPATH, _ = w4h.file_setup(well_data_path=well_data, verbose=verbose, log=log, **file_setup_kwargs)             
                headerDataPATH = well_metadata
            elif well_metadata is None:
                downholeDataPATH, _ = w4h.file_setup(well_data_path=well_data, verbose=verbose, log=log, **file_setup_kwargs)             

    elif isinstance(well_data, pd.DataFrame):
        if isinstance(well_metadata, pd.DataFrame):
            downholeDataPATH = well_data
            headerDataPATH = well_metadata
        elif isinstance(well_metadata, pathlib.PurePath) or isinstance(well_metadata, str):
            _, headerDataPATH = w4h.file_setup(well_data_path=well_metadata, metadata_path=well_metadata, verbose=verbose, log=log, **file_setup_kwargs)                
            downholeDataPATH = well_data
        else:
            print('ERROR: well_metadata must be a string filepath, a pathlib.Path object, or pandas.DataFrame')
    else:
        print('ERROR: well_data must be a string filepath, a pathlib.Path object, or pandas.DataFrame')

    #Get pandas dataframes from input
    read_raw_txt_kwargs = {k: v for k, v in locals().items() if k in w4h.read_raw_txt.__code__.co_varnames}
    downholeDataIN, headerDataIN = w4h.read_raw_txt(data_filepath=downholeDataPATH, metadata_filepath=headerDataPATH, verbose=verbose, log=log, **read_raw_txt_kwargs) #Functions to read data into dataframes. Also excludes extraneous columns, and drops header data with no location information

    #Define data types (file will need to be udpated)
    downholeData = w4h.define_dtypes(df=downholeDataIN, dtype_file='downholeDataTypes.txt', log=log)
    headerData = w4h.define_dtypes(df=headerDataIN, dtype_file='headerDataTypes.txt', log=log)

    #Get Study area
    read_study_area_kwargs = {k: v for k, v in locals().items() if k in w4h.read_study_area.__code__.co_varnames}
    if study_area is None:
        studyAreaIN = None
        use_study_area = False
    else:
        studyAreaIN = w4h.read_study_area(studyareapath=study_area, log=log, **read_study_area_kwargs)
        use_study_area = True

    #Get surfaces and grid(s)
    read_grid_kwargs = {k: v for k, v in locals().items() if k in w4h.read_grid_kwargs.__code__.co_varnames}

    modelGridPath = model_grid
    surfaceElevPath = surf_elev_file
    bedrockElevPath = bedrock_elev_file
    #UPDATE: allow other types of model grid read ***
    modelGrid = w4h.read_grid(datapath=modelGridPath, grid_type='model', study_area=studyAreaIN,  clip_to_studyarea=use_study_area, log=log, **read_grid_kwargs)
    surfaceElevGridIN = w4h.read_grid(datapath=surfaceElevPath, grid_type='surface', study_area=studyAreaIN, clip_to_studyarea=use_study_area, log=log, **read_grid_kwargs)
    bedrockElevGridIN = w4h.read_grid(datapath=bedrockElevPath, grid_type='bedrock', study_area=studyAreaIN, clip_to_studyarea=use_study_area, log=log, **read_grid_kwargs)

    #Add control points
    #UPDATE: Code here for adding in control points ***

    #UPDATE: MAKE SURE CRS's all align ***

    #Convert headerData to have geometry
    coords2geometry_kwargs = {k: v for k, v in locals().items() if k in w4h.coords2geometry.__code__.co_varnames}
    headerData = w4h.coords2geometry(df=headerData, xcol=xcol, ycol=ycol, zcol=zcol, log=log, **coords2geometry_kwargs)
    clip_gdf2study_area_kwargs = {k: v for k, v in locals().items() if k in w4h.clip_gdf2study_area.__code__.co_varnames}
    headerData = w4h.clip_gdf2study_area(study_area=studyAreaIN, gdf=headerData, log=log, **clip_gdf2study_area_kwargs)

    #Clean up data
    downholeData = w4h.remove_nonlocated(downholeData, headerData, log=log, verbose=verbose)
    headerData = w4h.remove_no_topo(df=headerData, zcol=zcol, verbose=verbose, log=log)

    drop_no_depth_kwargs = {k: v for k, v in locals().items() if k in w4h.drop_no_depth.__code__.co_varnames}
    donwholeData = w4h.drop_no_depth(downholeData, verbose=verbose, top_col=top_col, bottom_col=bottom_col, log=log, **drop_no_depth_kwargs) #Drop records with no depth information

    drop_bad_depth_kwargs = {k: v for k, v in locals().items() if k in w4h.drop_bad_depth.__code__.co_varnames}
    donwholeData = w4h.drop_bad_depth(downholeData, verbose=verbose, top_col=top_col, bottom_col=bottom_col, depth_type=depth_type, log=log, **drop_bad_depth_kwargs)#Drop records with bad depth information (i.e., top depth > bottom depth) (Also calculates thickness of each record)

    drop_no_formation_kwargs = {k: v for k, v in locals().items() if k in w4h.drop_no_formation.__code__.co_varnames}
    downholeData = w4h.drop_no_formation(downholeData, description_col=description_col, verbose=verbose, log=log, **drop_no_formation_kwargs)

    #CLASSIFICATION
    #Read dictionary definitions and classify
    #UPDATE: START HERE AGAIN, double checck and get kwargs for all functions ***
    get_search_terms_kwargs = {k: v for k, v in locals().items() if k in w4h.get_search_terms.__code__.co_varnames}
    specTermsPATH, startTermsPATH, wildcardTermsPATH, = w4h.get_search_terms(spec_dir=lith_dict, start_dir=lith_dict_start, wildcard_dir=lith_dict_wildcard, log=log, **get_search_terms_kwargs)
    read_dictionary_terms_kwargs = {k: v for k, v in locals().items() if k in w4h.read_dictionary_terms.__code__.co_varnames}
    specTerms = w4h.read_dictionary_terms(dict_file=specTermsPATH, log=log, **read_dictionary_terms_kwargs)
    startTerms = w4h.read_dictionary_terms(dict_file=startTermsPATH, log=log, **read_dictionary_terms_kwargs)
    wildcardTerms = w4h.read_dictionary_terms(dict_file=wildcardTermsPATH, log=log, **read_dictionary_terms_kwargs)

    #Clean up dictionary terms
    specTerms.drop_duplicates(subset='FORMATION', inplace=True)
    specTerms.reset_index(inplace=True, drop=True)

    startTerms.drop_duplicates(subset='FORMATION', inplace=True)
    startTerms.reset_index(inplace=True, drop=True)

    wildcardTerms.drop_duplicates(subset='FORMATION', inplace=True)
    wildcardTerms.reset_index(inplace=True, drop=True)

    #CLASSIFICATIONS
    #Exact match classifications
    downholeData = w4h.specific_define(downholeData, specTerms, verbose=verbose, log=log)
    
    #.startswith classifications
    classifedDF, searchDF = w4h.split_defined(downholeData)
    searchDF = w4h.start_define(df=searchDF, terms_df=startTerms, verbose=verbose, log=log)
    downholeData = w4h.remerge_data(classifieddf=classifedDF, searchdf=searchDF) #UPDATE: Needed? ***    

    #wildcard/any substring match classifications    
    classifedDF, searchDF = w4h.split_defined(downholeData)
    searchDF = w4h.wildcard_define(df=searchDF, terms_df=wildcardTerms, verbose=verbose, log=log)
    downholeData = w4h.remerge_data(classifieddf=classifedDF, searchdf=searchDF) #UPDATE: Needed? ***    
    
    #Depth classification
    classifedDF, searchDF = w4h.split_defined(downholeData)
    searchDF = w4h.depth_define(searchDF, thresh=550, verbose=verbose, log=log)
    downholeData = w4h.remerge_data(classifieddf=classifedDF, searchdf=searchDF) #UPDATE: Needed? ***
    
    #Fill unclassified data
    downholeData = w4h.fill_unclassified(downholeData)
    
    #Add target interpratations
    targetInterpDF = w4h.read_lithologies(log=log)
    downholeData = w4h.merge_lithologies(df=downholeData, targinterps_df=targetInterpDF)
    
    #Get ready for next steps
    wellsDF = w4h.get_unique_wells(downholeData)
    downholeData = w4h.sort_dataframe(df=downholeData, sort_cols=['API_NUMBER','TOP'], remove_nans=True)

    #Analyze Surface(s) and grid(s)
    bedrockGrid, surfaceGrid = w4h.align_rasters(grids_unaligned=[bedrockElevGridIN, surfaceElevGridIN], modelgrid=modelGrid, log=log)
    driftThickGrid, layerThickGrid = w4h.get_drift_thick(surface=surfaceGrid, bedrock=bedrockGrid, layers=9, plot=False, log=log)
    headerData = w4h.sample_raster_points(raster=bedrockGrid, points_df=headerData, new_col='BEDROCK_ELEV_FT', log=log)
    headerData = w4h.sample_raster_points(raster=surfaceGrid, points_df=headerData, new_col='SURFACE_ELEV_FT', log=log)
    headerData = w4h.sample_raster_points(raster=driftThickGrid, points_df=headerData, new_col='BEDROCK_DEPTH_FT', log=log)
    headerData = w4h.sample_raster_points(raster=layerThickGrid, points_df=headerData, new_col='LAYER_THICK_FT', log=log)
    headerData = w4h.get_layer_depths(well_metadata=headerData, no_layers=9, log=log)

    #UPDATE: Check if this actually works, I think they should be copies of each other if well_metadata is not specified and not found using the metadata_filename pattern in file_setup() ***
    #Merge header and data into one df, if applicable
    if downholeData.values.base is headerData.values.base:
        pass
    else:
        #downholeData = pd.merge(left = downholeData, right = headerData, on=idcol)
        downholeData = w4h.merge_tables(data_df=downholeData,  header_df=headerData, data_cols=None, header_cols=None, on=idcol, how='inner', auto_pick_cols=True, log=log)
    
    #downholeData = downholeData.copy()
    #UPDATE: Potentially need to remove duplicate columns here, I think I fixed that tho ***
    resdf = w4h.layer_target_thick(downholeData, layers=9, outfile_prefix='CoarseFine', log=log)
    layers_data = w4h.layer_interp(points=resdf, layers=9, grid=modelGrid, method='lin', log=log)

    if export_dir is None:
        if well_data.is_dir():
            export_dir = well_data.joinpath('Output')
        else:
            export_dir = well_data.parent.joinpath()
        
        if not export_dir.exists():
            try:
                export_dir.mkdir()
            except:
                pass

    w4h.export_grids(grid_data=layers_data, out_path=export_dir, file_id='',filetype='tif', variable_sep=True, date_stamp=True, log=log)
    return resdf, layers_data