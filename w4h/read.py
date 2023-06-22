import datetime
import inspect
import json
import os
import pathlib

import pandas as pd
import numpy as np

repoDir = pathlib.Path(os.getcwd())

from w4h import logger_function

# Gets the current date for use with in code
def get_current_date():
    """ Gets the current date to help with finding the most recent file
        ---------------------
        Parameters:
            None

        ---------------------
        Returns:
            todayDate   : datetime object with today's date
            dateSuffix  : str to use for naming output files
    """
    todayDate = datetime.date.today()
    todayDateStr = str(todayDate)
    dateSuffix = '_'+todayDateStr
    return todayDate, dateSuffix

#Function to get most recent file 
def get_most_recent(dir=str(repoDir)+'/resources', glob_pattern='*', verbose=True):
    """Function to find the most recent file with the indicated pattern, using pathlib.glob function.

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
    """
    todayDate = datetime.date.today()
    todayDateStr = str(todayDate)

    files = pathlib.Path(dir).rglob(glob_pattern) #Get all the files that fit the pattern
    fileDates = []
    for f in files: #Get the file dates from their file modification times
        fileDates.append(np.datetime64(datetime.datetime.fromtimestamp(os.path.getmtime(f))))
    
    if fileDates == []:
        #If no files found that match pattern, return an empty pathlib.Path()
        if verbose:
            print('No file found in {} matching {} pattern'.format(dir, glob_pattern))
        mostRecentFile = pathlib.Path()
        return mostRecentFile
    else:
        globInd = np.argmin(np.datetime64(todayDateStr) - np.array(fileDates)) #Find the index of the most recent file

    #Iterate through glob/files again (need to recreate glob)
    files = pathlib.Path(dir).rglob(glob_pattern)
    for j, f in enumerate(files):
        if j == globInd:
            mostRecentFile=f
            break
    
    if verbose:
        print('Most Recent version of this file is : '+mostRecentFile.name)

    return mostRecentFile

#Function to setup files of interest
def file_setup(well_data, metadata_path=None, data_filename='*ISGS_DOWNHOLE_DATA*.txt', metadata_filename='*ISGS_HEADER*.txt', log_dir=None, verbose=False, log=False):
    """Function to setup files, assuming data, metadata, and elevation/location are in separate files (there should be one "key"/identifying column consistent across all files to join/merge them later)

    This function may not be useful if files are organized differently than this structure. 
    If that is the case, it is recommended to use the get_most_recent() function for each individual file if needed.
    It may also be of use to simply skip this function altogether and directly define each filepath in a manner that can be used by pandas.read_csv()

    Parameters
    ----------
    well_data : str or pathlib.Path object
        Str or pathlib.Path to directory containing input files, by default str(repoDir)+'/resources'
    metadata_path : str or pathlib.Path object, optional
        Str or pathlib.Path to directory containing input metadata files, by default str(repoDir)+'/resources'
    data_filename : str, optional
        Pattern used by pathlib.glob() to get the most recent data file, by default '*ISGS_DOWNHOLE_DATA*.txt'
    metadata_filename : str, optional
        Pattern used by pathlib.glob() to get the most recent metadata file, by default '*ISGS_HEADER*.txt'
    verbose : bool, default = False
        Whether to print name of files to terminal, by default True
    log : bool, default = True
        Whether to log inputs and outputs to log file.

    Returns
    -------
    tuple
        Tuple with (well_data, metadata)
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    #Define  filepath variables to be used later for reading/writing files
    data_path = pathlib.Path(well_data)
    if metadata_path is None:
        origMetaPath = None
        metadata_path=data_path
    else:
        origMetaPath = metadata_path
        metadata_path=pathlib.Path(metadata_path)

    #If input path is a directory, find most recent version of the file. If file, just read the file
    if data_path.is_dir():
        downholeDataFILE = get_most_recent(data_path, data_filename, verbose=verbose)
    else:
        downholeDataFILE = data_path
    
    if metadata_path.is_dir():
        headerDataFILE = get_most_recent(metadata_path, metadata_filename, verbose=verbose)
        if headerDataFILE == []:
            headerDataFILE = downholeDataFILE
    else:
        if origMetaPath is None:
            headerDataFILE = downholeDataFILE
        else:
            headerDataFILE = metadata_path
       #Set all input as pathlib.Path objects (may be redundant, but just in case)
    downholeDataPATH = pathlib.Path(downholeDataFILE)
    headerDataPATH = pathlib.Path(headerDataFILE)

    if verbose:
        print('Using the following files:')
        print('\t', downholeDataFILE)
        print('\t', headerDataFILE)

    #Define datatypes, to use later
    #downholeDataDTYPES = {'ID':np.uint32, "API_NUMBER":np.uint64,"TABLE_NAME":str,"WHO":str,"INTERPRET_DATE":str,"FORMATION":str,"THICKNESS":np.float64,"TOP":np.float64,"BOTTOM":np.float64}
    #headerDataDTYPES = {'ID':np.uint32,'API_NUMBER':np.uint64,"TDFORMATION":str,"PRODFORM":str,"TOTAL_DEPTH":np.float64,"SECTION":np.float64,"TWP":np.float64,"TDIR":str,"RNG":np.float64,"RDIR":str,"MERIDIAN":np.float64,"FARM_NAME":str,"NSFOOT":np.float64,"NSDIR":str,"EWFOOT":np.float64,"EWDIR":str,"QUARTERS":str,"ELEVATION":np.float64,"ELEVREF":str,"COMP_DATE":str,"STATUS":str,"FARM_NUM":str,"COUNTY_CODE":np.float64,"PERMIT_NUMBER":str,"COMPANY_NAME":str,"COMPANY_CODE":str,"PERMIT_DATE":str,"CORNER":str,"LATITUDE":np.float64,"LONGITUDE":np.float64,"ENTERED_BY":str,"UPDDATE":str,"ELEVSOURCE":str, "ELEV_FT":np.float64}
    return downholeDataPATH, headerDataPATH

#Read raw data by text
def read_raw_txt(data_filepath, metadata_filepath, data_cols=None, metadata_cols=None, xcol='LONGITUDE', ycol='LATITUDE', id_col='API_NUMBER', encoding='latin-1', verbose=False, log=False):
    """Easy function to read raw .txt files output from (for example), an Access database

    Parameters
    ----------
    data_filepath : str
        Filename of the file containing data, including the extension.
    metadata_filepath : str
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
    log : bool, default = False
        Whether to log inputs and outputs to log file.

    Returns
    -------
    (pandas.DataFrame, pandas.DataFrame)
        Tuple/list with two pandas dataframes: (data, metadata)
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if data_cols is None:
        data_useCols = ["API_NUMBER","TABLE_NAME","FORMATION","THICKNESS","TOP","BOTTOM"]
    else:
        data_useCols= data_cols

    if metadata_cols is None:
        metadata_useCols = ['API_NUMBER',"TOTAL_DEPTH","SECTION","TWP","TDIR","RNG","RDIR","MERIDIAN","QUARTERS","ELEVATION","ELEVREF","COUNTY_CODE","LATITUDE","LONGITUDE","ELEVSOURCE"]
    else:
        metadata_useCols= metadata_cols

    #Check if input data is already 
    if isinstance(data_filepath, pd.DataFrame):
        downholeDataIN = data_filepath[data_useCols]
    else:
        downholeDataIN = pd.read_csv(data_filepath, sep=',', header='infer', encoding=encoding, usecols=data_useCols)
    
    if isinstance(metadata_filepath, pd.DataFrame):
        headerDataIN = metadata_filepath[metadata_useCols]
    else:
        headerDataIN = pd.read_csv(metadata_filepath, sep=',', header='infer', encoding=encoding, usecols=metadata_useCols)

    #Drop data with no API        
    downholeDataIN = downholeDataIN.dropna(subset=[id_col]) #Drop data with no API
    headerDataIN = headerDataIN.dropna(subset=[id_col]) #Drop metadata with no API

    #Drop data with no or missing location information
    headerDataIN = headerDataIN.dropna(subset=[ycol]) 
    headerDataIN = headerDataIN.dropna(subset=[xcol])
    
    #Reset index so index goes from 0 in numerical/integer order
    headerDataIN.reset_index(inplace=True, drop=True)
    downholeDataIN.reset_index(inplace=True, drop=True)
    
    #Print outputs to terminal, if designated
    if verbose:
        print('Data file has ' + str(downholeDataIN.shape[0])+' valid well records.')
        print("Metadata file has "+str(headerDataIN.shape[0])+" unique wells with valid location information.")
    
    return downholeDataIN, headerDataIN

#Read file with xyz data
def read_xyz(xyzpath, dtypes=None, verbose=False, log=False):
    """Function to read file containing xyz data (elevation/location)

    Parameters
    ----------
    xyzpath : str or pathlib.Path
        Filepath of the xyz file, including extension
    dtypes : dict, default = None
        Dictionary containing the datatypes for the columns int he xyz file. If None, {'ID':np.uint32,'API_NUMBER':np.uint64,'LATITUDE':np.float64,'LONGITUDE':np.float64,'ELEV_FT':np.float64}, by default None
    verbose : bool, default = False
        Whether to print the number of xyz records to the terminal, by default False
    log : bool, default = False
        Whether to log inputs and outputs to log file.

    Returns
    -------
    pandas.DataFrame
        Pandas dataframe containing the elevation and location data
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if dtypes is None:
        xyzDTypes = {'ID':np.uint32,'API_NUMBER':np.uint64,'LATITUDE':np.float64,'LONGITUDE':np.float64,'ELEV_FT':np.float64}

    xyzDataIN = pd.read_csv(xyzpath, sep=',', header='infer', dtype=xyzDTypes, index_col='ID')
    
    if verbose:
        print('XYZ file has ' + str(xyzDataIN.shape[0])+' records with elevation and location.')
    return xyzDataIN

#Get filepath of resource in resource folder
def __get_resource_path(res):
    """Function to get the path to the resource folder (used in other functions)

    Parameters
    ----------
    res : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    repoPath = pathlib.Path(__file__).parent.parent
    repoPathStr = str(repoPath).replace('\\', '/').replace('\\'[0], '/')
    resource = repoPathStr+'/resources/'+res
    return resource

#Read dictionary file into dictionary variable
def read_dict(file, keytype='np'):
    """Function to read a text file with a dictionary in it into a python dictionary

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
    """

    with open(file, 'r') as f:
        data= f.read()

    jsDict = json.loads(data)
    if keytype=='np':
        for k in jsDict.keys():
            jsDict[k] = getattr(np, jsDict[k])
    
    return jsDict

#Define the datatypes for a dataframe
def define_dtypes(df, dtypes=None, dtype_file=None, dtype_dir=str(repoDir)+'/resources/', log=False):
    """Function to define datatypes of a dataframe, especially with file-indicated dyptes

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
    log : bool, default = False
        Whether to log inputs and outputs to log file.

    Returns
    -------
    dfout : pandas.DataFrame
        Pandas dataframe containing redefined columns
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    dfout = df.copy()
    
    if dtypes is None:
        if isinstance(dtype_dir, pathlib.PurePath):
            dtype_dir = dtype_dir.as_posix()
        if dtype_dir[-1] != '/':
            dtype_dir = dtype_dir + '/'

        if dtype_file is None:
            print('ERROR: Either dtype_file (and dtype_dir) or dtypes must be defined')
            return 
        
        dtype_file = pathlib.Path(dtype_dir).joinpath(dtype_file)
        dtypes = read_dict(file=dtype_file)
    
    dfcols = dfout.columns
    for i in range(0, np.shape(dfout)[1]):
        dfout.iloc[:,i] = df.iloc[:,i].astype(dtypes[dfcols[i]])

    return dfout

#Define the search term filepaths
def get_search_terms(spec_dir=str(repoDir)+'/resources/', spec_glob_pattern='*SearchTerms-Specific*', 
                     start_dir=None, start_glob_pattern = '*SearchTerms-Start*', 
                     wildcard_dir=None, wildcard_glob_pattern='*SearchTerms-Wildcard',
                     log=False):
    """Read in dictionary files for downhole data

    Parameters
    ----------
    spec_dir : str or pathlib.Path, optional
        Directory where the file containing the specific search terms is located, by default str(repoDir)+'/resources/'
    spec_glob_pattern : str, optional
        Search string used by pathlib.glob() to find the most recent file of interest, uses get_most_recent() function, by default '*SearchTerms-Specific*'
    start_dir : str or None, optional
        Directory where the file containing the start search terms is located, by default None
    start_glob_pattern : str, optional
        Search string used by pathlib.glob() to find the most recent file of interest, uses get_most_recent() function, by default '*SearchTerms-Start*'
    wildcard_dir : str or pathlib.Path, default = None
        Directory where the file containing the wildcard search terms is located, by default None    
    wildcard_glob_pattern : str, default = '*SearchTerms-Wildcard'
        Search string used by pathlib.glob() to find the most recent file of interest, uses get_most_recent() function, by default '*SearchTerms-Wildcard*'
    log : bool, default = True
        Whether to log inputs and outputs to log file.        

    Returns
    -------
    (specTermsPath, startTermsPath, wilcardTermsPath) : tuple
        Tuple containing the pandas dataframes with specific search terms,  with start search terms, and with wildcard search terms
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)
    
    #specTermsFile = "SearchTerms-Specific_BedrockOrNo_2022-09.csv" #Specific matches
    #startTermsFile = "SearchTerms-Start_BedrockOrNo.csv" #Wildcard matches for the start of the description
    
    if start_dir is None:
        start_dir = spec_dir

    if wildcard_dir is None:
        if start_dir is None:
            wildcard_dir = spec_dir        
        else:
            wildcard_dir = start_dir        


    specTermsPath = get_most_recent(spec_dir, spec_glob_pattern)
    startTermsPath = get_most_recent(start_dir, start_glob_pattern)
    wilcardTermsPath = get_most_recent(wildcard_dir, wildcard_glob_pattern)
    
    return specTermsPath, startTermsPath, wilcardTermsPath

#Read files into pandas dataframes
def read_dictionary_terms(dict_file, id_col='ID', search_col='FORMATION', definition_col='INTERPRETATION', class_flag_col='CLASS_FLAG', dictionary_type=None, class_flag=6, rem_extra_cols=True, log=False):
    """Function to read dictionary terms from file into pandas dataframe

    Parameters
    ----------
    dict_file : str or pathlib.Path object, or list of these
        File or list of files to be read
    search_col : str, default = 'FORMATION'
        Name of column containing search terms (geologic formations)
    definition_col : str, default = 'INTERPRETATION'
        Name of column containing interpretations of search terms (lithologies)
    dictionary_type : str or None, {None, 'exact', 'start', 'wildcard',}
        Indicator of which kind of dictionary terms to be read in: None, 'exact', 'start', or 'wildcard' by default None.
            - If None, uses name of file to try to determine. If it cannot, it will default to using the classification flag from class_flag
            - If 'exact', will be used to search for exact matches to geologic descriptions
            - If 'start', will be used as with the .startswith() string method to find inexact matches to geologic descriptions
            - If 'wildcard', will be used to find any matching substring for inexact geologic matches
    class_flag : int, default = 1
        Classification flag to be used if dictionary_type is None and cannot be otherwise determined, by default 1
    rem_extra_cols : bool, default = True
        Whether to remove the extra columns from the input file after it is read in as a pandas dataframe, by default True
    log : bool, default = False
        Whether to log inputs and outputs to log file.

    Returns
    -------
    dict_terms : pandas.DataFrame
        Pandas dataframe with formatting ready to be used in the classification steps of this package
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    #Read files into pandas dataframes
    dict_terms = []
    if type(dict_file) is list:
        for f in dict_file:
            dict_terms.append(pd.read_csv(f))
            if id_col in dict_terms.columns:
                dict_terms.set_index(id_col, drop=True, inplace=True)
    else:
        dict_file = pathlib.Path(dict_file)
        if dict_file.exists():
            dict_terms.append(pd.read_csv(dict_file))
            if id_col in dict_terms[-1].columns:
                dict_terms[-1].set_index(id_col, drop=True, inplace=True)
            dict_file = [dict_file]
        else:
            print('ERROR: dict_file ({}) does not exist.'.format(dict_file))
            return

    #Rename important columns
    searchTermList = ['searchterm', 'search', 'exact']
    startTermList = ['startterm', 'start', 'startswith']
    wildcardTermList = ['wildcard', 'substring', ]

    #Recast all columns to datatypes of headerData to defined types
    dict_termDtypes = {search_col:str, definition_col:str, class_flag_col:np.uint8}

    if dictionary_type is None:
        dictionary_type='' #Allow string methods on this variable

    #Iterating, to allow reading of multiple dict file at once (also works with just one at at time)
    for i, d in enumerate(dict_terms):
        if dictionary_type.lower() in searchTermList or (dictionary_type=='' and 'spec' in str(dict_file[i]).lower()):
            d[class_flag_col] = 1
        elif dictionary_type.lower() in startTermList or (dictionary_type=='' and 'start' in str(dict_file[i]).lower()):
            d[class_flag_col] = 4 #Start term classification flag
        elif dictionary_type.lower() in wildcardTermList or (dictionary_type=='' and 'wildcard' in str(dict_file[i]).lower()):
            d[class_flag_col] = 5 #Wildcard term classification flag
        else:
            d[class_flag_col] = class_flag #Custom classification flag, defined as argument
        #1: exact classification match, 2: (not defined...ML?), 3: bedrock classification for obvious bedrock, 4: start term, 5: wildcard/substring, 6: Undefined

        #Rename columns so it is consistent through rest of code
        if search_col != 'FORMATION':
            d.rename(columns={search_col:'FORMATION'}, inplace=True)
        if definition_col != 'INTERPRETATION':
            d.rename(columns={definition_col:'INTERPRETATION'}, inplace=True)
        if class_flag_col != 'CLASS_FLAG':
            d.rename(columns={class_flag_col:'CLASS_FLAG'}, inplace=True)

        #Cast all columns as type str, if not already
        for i in range(0, np.shape(d)[1]):
            if d.iloc[:,i].name in list(dict_termDtypes.keys()):
                d.iloc[:,i] = d.iloc[:,i].astype(dict_termDtypes[d.iloc[:,i].name])
    
        #Delete duplicate definitions
        d.drop_duplicates(subset=search_col,inplace=True) #Apparently, there are some duplicate definitions, which need to be deleted first
        d.reset_index(inplace=True, drop=True)

    #If only one file designated, convert it so it's no longer a list, but just a dataframe
    if len(dict_terms) == 1:
        dict_terms = dict_terms[0]
    
    #Whether to remove extra columns that aren't needed from dataframe
    if rem_extra_cols:
        dict_terms = dict_terms[[search_col, definition_col, class_flag_col]]

    return dict_terms

#Function to read lithology file into pandas dataframe
def read_lithologies(litho_dir=None, lith_file=None, interp_col='LITHOLOGY', target_col='CODE', use_cols=None, log=False):
    """Function to read lithology file into pandas dataframe

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
    log : bool, default = True
        Whether to log inputs and outputs to log file.

    Returns
    -------
    pandas.DataFrame
        Pandas dataframe with lithology information
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    #dictDir = "\\\\isgs-sinkhole\\geophysics\\Balikian\\ISWS_HydroGeo\\WellDataAutoClassification\\SupportingDocs\\"
    if litho_dir is None:
        litho_dir=str(repoDir)+'/resources/'
    elif isinstance(litho_dir, pathlib.PurePath):
        litho_dir = litho_dir.as_posix()
    
    litho_dir.replace('\\', '/')
    litho_dir.replace('\\'[-1], '/')
    if litho_dir[-1] != '/':
        litho_dir = litho_dir+'/'

    if lith_file is None:
        lith_file='Lithology_Interp_FineCoarse.csv'
    
    if use_cols is None:
        use_cols = ['LITHOLOGY', 'CODE']

    lithFPath = pathlib.Path(litho_dir+lith_file)
    lithoDF = pd.read_csv(lithFPath, usecols=use_cols)

    lithoDF.rename(columns={interp_col:'INTERPRETATION', target_col:'TARGET'}, inplace=True)

    return lithoDF