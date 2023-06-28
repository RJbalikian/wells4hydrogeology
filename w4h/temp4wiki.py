import w4h
# Workflow

#Input Data
## file_setup()    
    w4h.file_setup(
       *well_data: str or pathlib.PurePath(),
       *metadata: str, pathlib.PurePath(), or None = None,
        data_filename: str = '*ISGS_DOWNHOLE_DATA*.txt',
        metadata_filename: str = '*ISGS_HEADER*.txt',
        log_dir: str, pathlib.PurePath(), or None = None,
        )

## read_raw_txt()
    w4h.read_raw_txt(
       +data_filepath, 
       +metadata_filepath, 
        data_cols=None, 
        metadata_cols=None, 
        xcol='LONGITUDE', 
        ycol='LATITUDE', 
        id_col='API_NUMBER', 
        encoding='latin-1'
        )

## define_dtypes()
    w4h.define_dtypes(
      +*undefined_df: pd.DataFrame,
        datatypes: dict, str, pathlib.PurePath, or None = None
        )
NOTE: define_dtypes() is repeated for both the well_data and well_metadata file, so this should be one file or dictionary containing datatypes for both datasets

## read_study_area()
    w4h.read_study_area(
        study_area_path: str or pathlib.Path
        study_area_crs: str = 'EPSG:4269',
        )

## read_grid()
    w4h.read_grid(
       *grid_path: None,
        grid_type: str = 'model',
        no_data_val: int = 0,
        use_service: {False, 'wms', 'wcs'}
       +study_area: Any | None = None,
        study_area_crs: Any | None = None,
        grid_crs: Any | None = None,
        **kwargs: Any
        )
NOTE: read_grid() is used for the model grid, surface elevation, and bedrock elevation

## coords2geometry()
    w4h.coords2geometry(    
       +df_no_geometry: pandas.DataFrame = metadata or well_data,
        xcol: str = 'LONGITUDE',
        ycol: str = 'LATITUDE',
        zcol: str = 'ELEV_FT',
        input_coords_crs: str = 'EPSG:4269',
        use_z: bool = False
        )
## clip_gdf2study_area()
    w4h.clip_gdf2study_area

# Clean Data
## remove_nonlocated()
    w4h.remove_nonlocated

## remove_no_topo()
    w4h.remove_no_topo

## drop_no_depth()
    w4h.drop_no_depth

## drop_bad_depth()
    w4h.drop_bad_depth

## drop_no_formation()
    w4h.drop_no_formation

# Classify Data
## get_search_terms()
    w4h.get_search_terms

## read_dictionary_terms()
    w4h.read_dictionary_terms
NOTE: w4h.read_dictionary_terms() is run three times:
- once for exact terms,
- once for starting terms, and
- once for wildcard terms

## specific_define()
    w4h.specific_define()

## start_define()
    w4h.split_defined
    w4h.start_define
    w4h.remerge_data

## wildcard_define()
    w4h.split_defined
    w4h.wildcard_define
    w4h.remerge_data

## depth_define()
    w4h.split_defined
    w4h.depth_define
    w4h.remerge_data

## fill_unclassified()
    w4h.fill_unclassified

## read_lithologies()
    w4h.read_lithologies

## merge_lithologies()
    w4h.merge_lithologies

# Raster Analysis
## align_raster()
    w4h.align_rasters

## get_drift_thick()
    w4h.get_drift_thick

## sample_raster_points()
    w4h.sample_raster_points
NOTE: sample_raster_points() is run 4 times to sample:
- The bedrock grid
- The surface grid
- The derived drift thickness grid
- The derived layer thigkness grid


# Model Layer Calculations
## get_layer_depths()
    w4h.get_layer_depths

## layer_target_thick()
    w4h.layer_target_thick

## layer_interp()
    w4h.layer_interp

## export_grids()
    w4h.export_grids

NOTE: the w4h.run() function returns a tuple (point_data, grid_data)