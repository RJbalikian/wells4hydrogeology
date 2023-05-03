Module w4h.read
===============

Functions
---------

    
`define_dtypes(df, dtypes=None, dtype_file=None, dtype_dir='C:\\Users\\riley\\LocalData\\Code\\Github\\wells4hydrogeology/resources/')`
:   Function to define datatypes of a dataframe, especially with file-indicated dyptes
    
    Parameters
    ----------
    df : pandas.DataFrame
        Pandas dataframe with columns whose datatypes need to be (re)defined
    dtypes : dict or None, default = None
        Dictionary containing datatypes, to be used in pandas.DataFrame.astype() function. If None, will read from file indicated by dtype_file (which must be defined, along with dtype_dir), by default None
    dtype_file : str or None, default = None
        Filename of file containing datatypes (txt file in dictionary format). If None, dtypes must be defined, by default None.
    dtype_dir : str or pathlib.Path obejct, default = str(repoDir)+'/resources/'
        Directory containing dtype_file, by default 
    
    Returns
    -------
    dfout : pandas.DataFrame
        Pandas dataframe containing redefined columns

    
`file_setup(db_dir, metadata_dir=None, xyz_dir=None, data_pattern='*ISGS_DOWNHOLE_DATA*.txt', metadata_pattern='*ISGS_HEADER*.txt', xyz_pattern='*xyzData*', verbose=False)`
:   Function to setup files, assuming data, metadata, and elevation/location are in separate files (there should be one "key"/identifying column consistent across all files to join/merge them later)
    
    This function may not be useful if files are organized differently than this structure. If that is the case, it is recommended to use the get_most_recent() function for each individual file needed.
    
    Parameters
    ----------
    db_dir : str or pathlib.Path object, optional
        Str or pathlib.Path to directory containing input files, by default str(repoDir)+'/resources'
    metadata_dir : str or pathlib.Path object, optional
        Str or pathlib.Path to directory containing input metadata files, by default str(repoDir)+'/resources'
    xyz_dir : str or pathlib.Path object, optional
        Str or pathlib.Path to directory containing input metadata files, by default str(repoDir)+'/resources'
    data_pattern : str, optional
        Pattern used by pathlib.glob() to get the most recent data file, by default '*ISGS_DOWNHOLE_DATA*.txt'
    metadata_pattern : str, optional
        Pattern used by pathlib.glob() to get the most recent metadata file, by default '*ISGS_HEADER*.txt'
    xyz_pattern : str, optional
        Pattern used by pathlib.glob() to get the most recent elevation/location file, by default '*xyzData*'
    verbose : bool, default = False
        Whether to print name of files to terminal, by default True
    
    Returns
    -------
    _type_
        _description_

    
`get_current_date()`
:   Gets the current date to help with finding the most recent file
    ---------------------
    Parameters:
        None
    
    ---------------------
    Returns:
        todayDate   : datetime object with today's date
        dateSuffix  : str to use for naming output files

    
`get_most_recent(dir='C:\\Users\\riley\\LocalData\\Code\\Github\\wells4hydrogeology/resources', glob_pattern='*', verbose=True)`
:   Function to find the most recent file with the indicated pattern, using pathlib.glob function.
    
    Parameters
    ----------
    dir : str or pathlib.Path object, optional
        Directory in which to find the most recent file, by default str(repoDir)+'/resources'
    glob_pattern : str, optional
        String used by the pathlib.glob() function/method for searching, by default '*'
    
    Returns
    -------
    pathlib.Path object
        Pathlib Path object of the most recent file fitting the glob pattern indicated in the glob_pattern parameter.

    
`get_search_terms(spec_dir='C:\\Users\\riley\\LocalData\\Code\\Github\\wells4hydrogeology/resources/', specStartPattern='*SearchTerms-Specific*', start_dir=None, startGlobPattern='*SearchTerms-Start*')`
:   Read in dictionary files for downhole data
    
    Parameters
    ----------
    spec_dir : str or pathlib.Path, optional
        Directory where the file containing the specific search terms is located, by default str(repoDir)+'/resources/'
    specStartPattern : str, optional
        Search string used by pathlib.glob() to find the most recent file of interest, uses get_most_recent() function, by default '*SearchTerms-Specific*'
    start_dir : str or None, optional
        Directory where the file containing the start search terms is located, by default None
    startGlobPattern : str, optional
        Search string used by pathlib.glob() to find the most recent file of interest, uses get_most_recent() function, by default '*SearchTerms-Start*'
    
    Returns
    -------
    (specTermsPath, startTermsPath) : tuple
        Tuple containing the pandas dataframe with specific search terms and one with start search terms

    
`read_dict(file, keytype='np')`
:   Function to read a text file with a dictionary in it into a python dictionary
    
    Parameters
    ----------
    file : str or pathlib.Path object
        Filepath to the file of interest containing the dictionary text
    keytype : str, optional
        String indicating the datatypes used in the text, currently only 'np' is implemented, by default 'np'
    
    Returns
    -------
    dict
        Dictionary translated from text file.

    
`read_dictionary_terms(dict_file, cols=None, col_types=None, dictionary_type=None, class_flag=1, rem_extra_cols=True)`
:   Function to read dictionary terms from file into pandas dataframe
    
    Parameters
    ----------
    dict_file : str or pathlib.Path object, or list of these
        File or list of files to be read
    cols : dict or None, default = None
        Dictionary containing columns to be renamed. If None, see source code for renaming actions, by default None
    col_types : dict or None, default = None
        Dictionary containing column types to be set. If None, see source code for renaming actions, by default None, by default None
    dictionary_type : str or None, {None, 'exact', 'start'}
        Indicator of which kind of dictionary terms to be read in: None, 'exact' or 'start', by default None.
            - If None, uses name of file to try to determine. If it cannot, it will default to using the classification flag from class_flag
            - If 'exact', will be used to search for exact matches to geologic descriptions
            - If 'start', will be used as with the .startswith() string method to find inexact matches to geologic descriptions
    class_flag : int, default = 1
        Classification flag to be used if dictionary_type is None and cannot be otherwise determined, by default 1
    rem_extra_cols : bool, default = True
        Whether to remove the extra columns from the input file after it is read in as a pandas dataframe, by default True
    
    Returns
    -------
    dict_terms : pandas.DataFrame
        Pandas dataframe with formatting ready to be used in the classification steps of this package

    
`read_lithologies(litho_dir=None, lith_file=None, interp_col='LITHOLOGY', target_col='CODE', use_cols=None)`
:   Function to read lithology file into pandas dataframe
    
    Parameters
    ----------
    litho_dir : str or pathlib.Path object, default = None
        Directory where lithology file is located. If None, default is in source coude, by default None
    lith_file : str, default = None
        Filename of lithology file. If None, default is in source coude, by default None
    interp_col : str, default = 'LITHOLOGY'
        Column to used to match interpretations
    target_col : str, default = 'CODE'
        Column to be used as target code
    use_cols : list, default = None
        Which columns to use when reading in dataframe. If None, defaults to ['LITHOLOGY', 'CODE'].
    Returns
    -------
    pandas.DataFrame
        Pandas dataframe with lithology information

    
`read_raw_txt(raw_dir, data_filename, metadata_filename, data_cols=None, metadata_cols=None, x_col='LONGITUDE', ycol='LATITUDE', id_col='API_NUMBER', encoding='latin-1', verbose=False)`
:   Easy function to read raw .txt files output from (for example), an Access database
    
    Parameters
    ----------
    raw_dir : str or pathlib.Path object
        String or pathlib.Path to directory containing the files.
    data_filename : str
        Filename of the file containing data, including the extension.
    metadata_filename : str
        Filename of the file containing metadata, including the extension.
    data_cols : list, default = None
        List with strings with names of columns from txt file to keep after reading. If None, ["API_NUMBER","TABLE_NAME","FORMATION","THICKNESS","TOP","BOTTOM"], by default None.
    metadata_cols : list, default = None
        List with strings with names of columns from txt file to keep after reading. If None, ['API_NUMBER',"TOTAL_DEPTH","SECTION","TWP","TDIR","RNG","RDIR","MERIDIAN","QUARTERS","ELEVATION","ELEVREF","COUNTY_CODE","LATITUDE","LONGITUDE","ELEVSOURCE"], by default None
    x_col : str, default = 'LONGITUDE'
        Name of column in metadata file indicating the x-location of the well, by default 'LONGITUDE'
    ycol : str, default = 'LATITUDE'
        Name of the column in metadata file indicating the y-location of the well, by default 'LATITUDE'
    id_col : str, default = 'API_NUMBER'
        Name of the column with the key/identifier that will be used to merge data later, by default 'API_NUMBER'
    encoding : str, default = 'latin-1'
        Encoding of the data in the input files, by default 'latin-1'
    verbose : bool, default = False
        Whether to print the number of rows in the input columns, by default False
    
    Returns
    -------
    (pandas.DataFrame, pandas.DataFrame)
        Tuple/list with two pandas dataframes: (data, metadata)

    
`read_xyz(rawdir, xyzfile, dtypes=None, verbose=False)`
:   Function to read file containing xyz data (elevation/location)
    
    Parameters
    ----------
    rawdir : str or pathlib.Path object
        String to the directory in which the xyz file is contained
    xyzfile : str
        String with the filename of the xyz file, including extension.
    dtypes : dict, default = None
        Dictionary containing the datatypes for the columns int he xyz file. If None, {'ID':np.uint32,'API_NUMBER':np.uint64,'LATITUDE':np.float64,'LONGITUDE':np.float64,'ELEV_FT':np.float64}, by default None
    verbose : bool, default = False
        Whether to print the number of xyz records to the terminal, by default False
    
    Returns
    -------
    pandas.DataFrame
        Pandas dataframe containing the elevation and location data