import w4h
# Workflow

## Read Data Inputs
#### file_setup()
file_setup() takes the inputs for the data files, checks the filepaths, and can get the most recent version of a file if multiple are located in a directory (if directory is specified)

    w4h.file_setup(
       *well_data: str or pathlib.PurePath(),
       *metadata: str, pathlib.PurePath(), or None = None,
        data_filename: str = '*ISGS_DOWNHOLE_DATA*.txt',
        metadata_filename: str = '*ISGS_HEADER*.txt',
        log_dir: str, pathlib.PurePath(), or None = None,
        )

#### read_raw_csv()
read_raw_csv() uses the filepaths specified as inputs and checked using file_setup() to read in the files as pandas DataFrames.

    w4h.read_raw_csv(
       +data_filepath, 
       +metadata_filepath, 
        data_cols=None, 
        metadata_cols=None, 
        xcol='LONGITUDE', 
        ycol='LATITUDE', 
        id_col='API_NUMBER', 
        encoding='latin-1',
        ...,
        **read_csv_kwargs = kwargs to be passed to pd.read_csv()
        )

#### define_dtypes()
define_dtypes() defines the datatypes of a dataframe, especially with file-indicated dtypes (but also with a dict specifying filetypes: {col:np.dtype})

    w4h.define_dtypes(
      +*undefined_df: pd.DataFrame,
        datatypes: dict, str, pathlib.PurePath, or None = None,
        ...
        )
    
> NOTE: define_dtypes() is repeated for both the well_data and well_metadata file, so this should be one file or dictionary containing datatypes for both datasets

#### read_study_area()
read_study_area() reads a study area geospatial file into a geopandas geodataframe.

    w4h.read_study_area(
        study_area_path: str or pathlib.Path
        study_area_crs: str = 'EPSG:4269',
        ...
        )

#### read_grid()
read_grid() reads a raster file into xarray using rioxarray/rasterio. For the model grid, other specifications can be made as well, for example, to manually create a grid.

    w4h.read_grid(
       *grid_path: None,
        grid_type: str = 'model',
        no_data_val: int = 0,
        use_service: {False, 'wms', 'wcs'}
       +study_area: Any | None = None,
        study_area_crs: Any | None = None,
        grid_crs: Any | None = None,
        ...,
        **kwargs: read_model_grid() kwargs
        )
> NOTE: read_grid() is used for the model grid, surface elevation, and bedrock elevation

#### coords2geometry()
coords2geometry() creates a geometry column for a geopandas.GeoDataFrame  in the specified coordinate reference system from a pandas.DataFrame with xy coordinates.

    w4h.coords2geometry(    
       +df_no_geometry: pandas.DataFrame = metadata or well_data,
        xcol: str = 'LONGITUDE',
        ycol: str = 'LATITUDE',
        zcol: str = 'ELEV_FT',
        input_coords_crs: str = 'EPSG:4269',
        use_z: bool = False
        ...
        )
#### clip_gdf2study_area()
clip_gdf2study_area() clips a geopandas.GeoDataFrame to only include features within study area.

    w4h.clip_gdf2study_area(
       +study_area: Any,
       +gdf: Any,
        input_coords_crs: str = 'EPSG:4269',
        ...
        )
    
## Clean Data
#### remove_nonlocated()
remove_nonlocated() is a function to remove wells and well intervals where there is no location information.

    w4h.remove_nonlocated(
       +df: Any,
       +metadata_df: Any,
       ...
        )
    
#### remove_no_topo()
remove_no_topo() removes wells that do not have topography data (topography is needed for layer selection later). This function is intended to be run on the metadata dataframe (or whichever dataframe has elevations) after elevations been added.

    w4h.remove_no_topo(
        +df: Any,
        zcol: str = 'ELEV_FT',
        no_data_val: str = '',
        ...
        )

#### remove_no_depth()
remove_no_depth() removes well intervals with no depth/elevation information.

    w4h.remove_no_depth(
       +df: Any,
        top_col: str = 'TOP',
        bottom_col: str = 'BOTTOM',
        no_data_val: str = '',
        ...
        )

#### remove_bad_depth()
remove_bad_depth() removes all records in the DataFrame with well interpretations where the depth information is bad (i.e., where the bottom of the record is neerer to the surface than the top)

    w4h.remove_bad_depth(
       +df: Any,
        top_col: str = 'TOP',
        bottom_col: str = 'BOTTOM',
        depth_type: str = 'depth',
        ...
        )
    
#### remove_no_formation()
remove_no_formation() removes all records in the DataFrame containing the well descriptions where no description is actually given.

    w4h.remove_no_formation(
       +df: Any,
        description_col: str = 'FORMATION',
        no_data_val: str = '',
        ...
        )

## Classify Data
#### get_search_terms()
get_search_terms() checks and gets the filepaths of the dictionary file(s) with classification terms for classifying well intervals later

    w4h.get_search_terms(
       *spec_path: str = str(repoDir) + '/resources/',
        spec_glob_pattern: str = '*SearchTerms-Specific*',
        start_path: Any | None = None,
        start_glob_pattern: str = '*SearchTerms-Start*',
        wildcard_path: Any | None = None,
        wildcard_glob_pattern: str = '*SearchTerms-Wildcard',
        ...
        )

#### read_dictionary_terms()
read_dictionary_terms() is a function to read dictionary terms from a filepath specified in get_search_terms() as a pandas.DataFrame

    w4h.read_dictionary_terms(
       +dict_file: Any,
        id_col: str = 'ID',
        search_col: str = 'DESCRIPTION',
        definition_col: str = 'LITHOLOGY',
        class_flag_col: str = 'CLASS_FLAG',
        dictionary_type: {None, 'exact', 'start', 'wildcard'},
        class_flag: int = 6,
        rem_extra_cols: bool = True,
        ...
        )
    
> NOTE: w4h.read_dictionary_terms() is run three times:
> - once for exact terms,
> - once for starting terms, and
> - once for wildcard terms

#### specific_define()
specific_define() classifies terms that have been specifically defined as exact matches in the terms_df DataFrame.

    w4h.specific_define(
        +df: Any,
        +terms_df: Any,
        description_col: str = 'FORMATION',
        terms_col: str = 'DESCRIPTION',
        ...
        )

#### start_define()
start_define() classifies descriptions where a substring of variable length that begins an interval description is matched against a predefined substring.
split_defined() and remerge_data() are used to split the main dataframe into those intervals that have been combined and those that have not (start_define() is performed on the latter), 
then re-merges the data back together using remerge_data()

    w4h.split_defined(+df: Any, classification_col: str = 'CLASS_FLAG', ...)
    w4h.start_define(
       +df: Any,
       +terms_df: Any,
        description_col: str = 'FORMATION',
        terms_col: str = 'DESCRIPTION',
        ...
        )
    w4h.remerge_data(+classifieddf: Any, +searchdf: Any)

#### wildcard_define()
wildcard_define() classifies descriptions where a substring of variable length *anywhere* in an interval description is matched against a predefined substring.
split_defined() and remerge_data() are used to split the main dataframe into those intervals that have been combined and those that have not (wildcard_define() is performed on the latter), 
then re-merges the data back together using remerge_data()

    w4h.split_defined(+df: Any, classification_col: str = 'CLASS_FLAG', ...)
    w4h.wildcard_define(
       +df: Any,
       +terms_df: Any,
        description_col: str = 'FORMATION',
        terms_col: str = 'DESCRIPTION',
        ...
        )
    w4h.remerge_data(+classifieddf: Any, +searchdf: Any)

#### depth_define()
depth_define() defines all intervals as bedrock that are lower than the depth or elevation defined by the parameter *thresh*
split_defined() and remerge_data() are used to split the main dataframe into those intervals that have been combined and those that have not (wildcard_define() is performed on the latter), 
then re-merges the data back together using remerge_data()

    w4h.split_defined(+df: Any, classification_col: str = 'CLASS_FLAG', ...)
    w4h.depth_define(
       +dfIN: Any,
        top_col: str = 'TOP',
        thresh: float = 550,
        ...
        )
    w4h.remerge_data(+classifieddf: Any, +searchdf: Any)

#### fill_unclassified()
fill_unclassified() fills empty rows in 'CLASS_FLAG' column (i.e., unclassified rows) with np.nan

    w4h.fill_unclassified(
       +df: Any,
        classification_col: str = 'CLASS_FLAG'
        )

#### read_lithologies()
read_lithologies() reads a lithology file into a pandas dataframe

    w4h.read_lithologies(
       *lith_file: Any | None = None,
        interp_col: str = 'LITHOLOGY',
        target_col: str = 'CODE',
        use_cols: Any | None = None,
        ...
        )

#### merge_lithologies()
merge_lithologies() merges lithologies and target booleans based on previous classifications

    w4h.merge_lithologies(
       +df: Any,
       +targinterps_df: Any,
        target_col: str = 'TARGET',
        target_class: str = 'bool'
        )

## Raster Analysis
#### align_raster()
align_rasters() reprojects two rasters and aligns their pixels with a model grid

    w4h.align_rasters(
       +grids_unaligned: Any,
       +modelgrid: Any,
        no_data_val: int = 0,
        ...
        )

#### get_drift_thick()
get_drift_thick() finds the distance from surface to bedrock and then divides that distance by the specified number of layers to geta model layer thickness.

    w4h.get_drift_thick(
       +surface: Any,
       +bedrock: Any,
       *layers: int = 9,
        plot: bool = False,
        ...
        )

#### sample_raster_points()
sample_raster_points() samples raster values in xarray to points from a geopandas geodataframe.

    w4h.sample_raster_points(
       +raster: Any,
       +points_df: Any,
        xcol: str = 'LONGITUDE',
        ycol: str = 'LATITUDE',
        new_col: str = 'SAMPLED',
        ...
        )
    
> NOTE: sample_raster_points() is run 4 times to sample the following grids at each well location:
> - The bedrock grid
> - The surface grid
> - The derived drift thickness grid
> - The derived layer thigkness grid

## Model Layer Calculations
#### get_layer_depths()
get_layer_depths() calculates depths and elevations of each model layer at each well based on surface elevation, bedrock elevation, and number of layers/layer thickness

    w4h.get_layer_depths(
       +well_metadata: Any,
        no_layers: int = 9,
        ...
        ) 

#### layer_target_thick()
layer_target_thick calculates the thickness of the target material in each layer at each well point. **This function does the core calculation at the heart of the entire process!**

    w4h.layer_target_thick(
       +df: Any,
        layers: int = 9,
        export_dir: Any | None = None,
        return_all: bool = False,
        outfile_prefix: str = '',
        depth_top_col: str = 'TOP',
        depth_bot_col: str = 'BOTTOM',
        ...
        )

#### layer_interp()
layer_interp() interpolates layer thickness results across the entire well field/study area, going from points to grid data. Uses the scipy.interpolate module.

    w4h.layer_interp(
       +points: Any,
       +grid: xr.DataArray or xr.DataSet,
        layers: Any | None = None,
        method: str = 'nearest',
        return_type: str = 'dataarray',
        export_dir: Any | None = None,
       +targetcol: str = 'TARG_THICK_PER',
       +lyrcol: str = 'LAYER',
        xcol: Any | None = None,
        ycol: Any | None = None,
        xcoord: str = 'x',
        ycoord: str = 'y',
        ...,
        **kwargs: scipy.interpolate kwargs determined by method argument
        )

#### export_grids()
export_grids() exports all grids (interpolated layer thickness, surface, bedrock, drift thickness, layer thickness) in specified format to specified location.
    w4h.export_grids(
       +grid_data: Any,
        out_path: Any,
        file_id: str = '',
        filetype: str = 'tif',
        variable_sep: bool = True,
        date_stamp: bool = True,
        ...
        )

> NOTE: the w4h.run() function returns a tuple (point_data, grid_data)