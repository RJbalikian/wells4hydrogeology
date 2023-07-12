import datetime
import inspect

import pandas as pd
import numpy as np

from w4h import logger_function
#The following flags are used to mark the classification method:
#- 0: Not classified
#- 1: Specific Search Term Match
#- 2: wPermits bedrock top pick
#- 3: Intervals >550' below ground surface
#- 4: Wildcard match (startTerm) - no context
#- 5: Wildcard match (any substring) - more liberal
#- Top of well?

#Define records with full search term

def specific_define(df, terms_df, description_col='FORMATION', terms_col='DESCRIPTION', verbose=False, log=False):
    """Function to classify terms that have been specifically defined in the terms_df.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe with unclassified well descriptions.
    terms_df : pandas.DataFrame
        Dataframe containing the classifications
    description_col : str, default='FORMATION'
        Column name in df containing the well descriptions, by default 'FORMATION'.
    terms_col : str, default='FORMATION'
        Column name in terms_df containing the classified descriptions, by default 'FORMATION'.
    verbose : bool, default=False
        Whether to print up results, by default False.

    Returns
    -------
    df_Interps : pandas.DataFrame
        Dataframe containing the well descriptions and their matched classifications.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if description_col != terms_col:
        terms_df.rename(columns={terms_col:description_col}, inplace=True)
        terms_col = description_col

    df[description_col] = df[description_col].astype(str)
    terms_df[terms_col] = terms_df[terms_col].astype(str)

    df[description_col] = df[description_col].str.casefold()
    terms_df[terms_col] = terms_df[terms_col].str.casefold()
    #df['FORMATION'] = df['FORMATION'].str.strip(['.,:?\t\s'])
    #terms_df['FORMATION'] = terms_df['FORMATION'].str.strip(['.,:?\t\s'])

    terms_df.drop_duplicates(subset=terms_col, keep='last', inplace=True)
    terms_df.reset_index(drop=True, inplace=True)
    
    df_Interps = pd.merge(df, terms_df.set_index(terms_col), on=description_col, how='left')
    df_Interps.rename(columns={description_col:'FORMATION'}, inplace=True)
    df_Interps['BEDROCK_FLAG'] = df_Interps['LITHOLOGY'] == 'BEDROCK'
    
    if verbose:
        print("Records Classified with full search term: "+str(int(df_Interps[df_Interps['CLASS_FLAG']==1]['CLASS_FLAG'].sum())))
        print("Records Classified with full search term: "+str(round((df_Interps[df_Interps['CLASS_FLAG']==1]['CLASS_FLAG'].sum()/df_Interps.shape[0])*100,2))+"% of data")

    return df_Interps


def split_defined(df, classification_col='CLASS_FLAG', verbose=False, log=False):
    """Function to split dataframe with well descriptions into two dataframes based on whether a row has been classified.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing all the well descriptions
    classification_col : str, default = 'CLASS_FLAG'
        Name of column containing the classification flag, by default 'CLASS_FLAG'
    verbose : bool, default = False
        Whether to print results, by default False
    log : bool, default = False
        Whether to log results to log file

    Returns
    -------
    Two-item tuple of pandas.Dataframe
        tuple[0] is dataframe containing classified data, tuple[1] is dataframe containing unclassified data.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    classifedDF= df[df[classification_col].notna()] #Already-classifed data
    searchDF = df[df[classification_col].isna()] #Unclassified data
    
    if verbose:
        print(str(searchDF.shape[0])+' records unclassified; isolated into searchDF.')
    
    return classifedDF, searchDF

#Classify downhole data by the initial substring
def start_define(df, terms_df, description_col='FORMATION', terms_col='DESCRIPTION', verbose=False, log=False):
    """Function to classify descriptions according to starting substring. 

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing all the well descriptions
    terms_df : pandas.DataFrame
        Dataframe containing all the startswith substrings to use for searching
    description_col : str, default = 'FORMATION'
        Name of column in df containing descriptions, by default 'FORMATION'
    terms_col : str, default = 'FORMATION'
        Name of column in terms_df containing startswith substring to match with description_col, by default 'FORMATION'
    verbose : bool, default = False
        Whether to print out results, by default False
    log : bool, default = True
        Whether to log results to log file

    Returns
    -------
    df : pandas.DataFrame
        Dataframe containing the original data and new classifications
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if verbose:
        #Estimate when it will end, based on test run
        estTime = df.shape[0]/3054409 * 6 #It took about 6 minutes to classify data with entire dataframe. This estimates the fraction of that it will take
        nowTime = datetime.datetime.now()
        endTime = nowTime+datetime.timedelta(minutes=estTime)
        print("Start Term process should be done by {:d}:{:02d}".format(endTime.hour, endTime.minute))

    #First, for each startterm, find all results in df that start with, add classification flag, and add interpretation.
    for i,s in enumerate(terms_df[terms_col]):
        df['CLASS_FLAG'].where(~df[description_col].str.startswith(s,na=False),4,inplace=True)
        df['LITHOLOGY'].where(~df[description_col].str.startswith(s,na=False),terms_df.loc[i,'LITHOLOGY'],inplace=True)
    df['BEDROCK_FLAG'].loc[df["LITHOLOGY"] == 'BEDROCK']
    
    if verbose:
        print("Records classified with start search term: "+str(int(df['CLASS_FLAG'].count())))
        print("Records classified with start search term: "+str(round((df['CLASS_FLAG'].count()/df.shape[0])*100,2))+"% of remaining data")
        #print("Records classified with both search terms: "+str(round(((df['CLASS_FLAG'].count()+specDF['CLASS_FLAG'].count())/downholeData_Interps.shape[0])*100,2))+"% of all data")
        #This step usually takes about 5-6 minutes
    return df

#Classify downhole data by any substring
def wildcard_define(df, terms_df, description_col='FORMATION', terms_col='DESCRIPTION', verbose=False, log=False):
    """Function to classify descriptions according to any substring. 

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing all the well descriptions
    terms_df : pandas.DataFrame
        Dataframe containing all the startswith substrings to use for searching
    description_col : str, default = 'FORMATION'
        Name of column in df containing descriptions, by default 'FORMATION'
    terms_col : str, default = 'FORMATION'
        Name of column in terms_df containing startswith substring to match with description_col, by default 'FORMATION'
    verbose : bool, default = False
        Whether to print out results, by default False
    log : bool, default = True
        Whether to log results to log file

    Returns
    -------
    df : pandas.DataFrame
        Dataframe containing the original data and new classifications
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if verbose:
        #Estimate when it will end, based on test run
        estTime = df.shape[0]/3054409 * 6 #It took about 6 minutes to classify data with entire dataframe. This estimates the fraction of that it will take
        nowTime = datetime.datetime.now()
        endTime = nowTime+datetime.timedelta(minutes=estTime)
        print("Wildcard Term process should be done by (?) {:d}:{:02d}".format(endTime.hour, endTime.minute))

    #First, for each startterm, find all results in df that start with, add classification flag, and add interpretation.
    for i,s in enumerate(terms_df[terms_col]):
        df['CLASS_FLAG'].where(~df[description_col].str.contains(s, case=False, regex=False, na=False), 5, inplace=True)
        df['LITHOLOGY'].where(~df[description_col].str.contains(s, case=False, regex=False, na=False),terms_df.loc[i,'LITHOLOGY'],inplace=True)
    df['BEDROCK_FLAG'].loc[df["LITHOLOGY"] == 'BEDROCK']
    
    if verbose:
        print("Records classified with wildcard search term: "+str(int(df['CLASS_FLAG'].count())))
        print("Records classified with wildcard search term: "+str(round((df['CLASS_FLAG'].count()/df.shape[0])*100,2))+"% of remaining data")
        #print("Records classified with both search terms: "+str(round(((df['CLASS_FLAG'].count()+specDF['CLASS_FLAG'].count())/downholeData_Interps.shape[0])*100,2))+"% of all data")
        #This step usually takes about 5-6 minutes
    return df

#Merge data back together
def remerge_data(classifieddf, searchdf):
    """Function to merge newly-classified (or not) and previously classified data

    Parameters
    ----------
    classifieddf : pandas.DataFrame
        Dataframe that had already been classified previously
    searchdf : pandas.DataFrame
        Dataframe with new classifications

    Returns
    -------
    remergeDF : pandas.DataFrame
        Dataframe containing all the data, merged back together
    """
    remergeDF = pd.concat([classifieddf,searchdf], join='inner').sort_index()
    return remergeDF

#Define well intervals by depth

def depth_define(dfIN, top_col='TOP', thresh=550.0, verbose=False, log=False):
    """Function to define all intervals lower than thresh as bedrock

    Parameters
    ----------
    dfIN : pandas.DataFrame
        Dataframe to classify
    top_col : str, default = 'TOP'
        Name of column that contains the depth information, likely of the top of the well interval, by default 'TOP'
    thresh : float, default = 550.0
        Depth (in units used in dfIN['top_col']) below which all intervals will be classified as bedrock, by default 550.0.
    verbose : bool, default = False
        Whether to print results, by default False
    log : bool, default = True
        Whether to log results to log file

    Returns
    -------
    df : pandas.DataFrame
        Dataframe containing intervals classified as bedrock due to depth
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    df = dfIN.copy()
    df['CLASS_FLAG'].mask(df[top_col]>thresh, 3 ,inplace=True) #Add a Classification Flag of 3 (bedrock b/c it's deepter than 550') to all records where the top of the interval is >550'
    df['BEDROCK_FLAG'].mask(df[top_col]>thresh, True, inplace=True)

    if verbose:
        if df.CLASS_FLAG.notnull().sum() == 0:
            brDepthClass = 0
        else:
            brDepthClass = df['CLASS_FLAG'].value_counts()[3.0]
        total = dfIN.shape[0]
        print("Records classified as bedrock that were deeper than "+str(thresh)+ "': " + str(brDepthClass))
        print("This represents "+str(round((brDepthClass)*100/total,2))+"% of the unclassified data in this dataframe.")
        
    return df

#Output data that still needs to be defined
def export_undefined(df, outdir):
    """Function to export terms that still need to be defined.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing at least some unclassified data
    outdir : str or pathlib.Path
        Directory to save file. Filename will be generated automatically based on today's date.

    Returns
    -------
    stillNeededDF : pandas.DataFrame
        Dataframe containing only unclassified terms, and the number of times they occur
    """
    import pathlib
    if isinstance(outdir, pathlib.PurePath):
        if not outdir.is_dir() or not outdir.exists():
            print('Please specify a valid directory for export. Filename is generated automatically.')
            return
        outdir = outdir.as_posix()
    else:
        outdir.replace('\\','/')
        outdir.replace('\\'[-1], '/')

    #Get directory path correct        
    if outdir[-1] != '/':
        outdir = outdir+'/'

    todayDate = datetime.date.today()
    todayDateStr = str(todayDate)
    searchDF = df[df['CLASS_FLAG'].isna()]
    
    stillNeededDF=searchDF['FORMATION'].value_counts()
    stillNeededDF.to_csv(outdir+'Undefined_'+todayDateStr+'.csv')
    return stillNeededDF

#Fill in unclassified rows' flags with 0
def fill_unclassified(df, classification_col='CLASS_FLAG'):
    """Fills unclassified rows in 'CLASS_FLAG' column with np.nan

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe on which to perform operation

    Returns
    -------
    df : pandas.DataFrame
        Dataframe on which operation has been performed
    """
    df[classification_col].fillna(0, inplace=True)
    return df

#Merge lithologies to main df based on classifications
def merge_lithologies(df, targinterps_df, interp_col='INTERPRETATION', target_col='TARGET', target_class='bool'):
    """Function to merge lithologies and target booleans based on classifications
    
    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing classified well data
    targinterps_df : pandas.DataFrame
        Dataframe containing lithologies and their target interpretations, depending on what the target is for this analysis (often, coarse materials=1, fine=0)
    target_col : str, default = 'TARGET'
        Name of column in targinterps_df containing the target interpretations
    target_class, default = 'bool'
        Whether the input column is using boolean values as its target indicator
        
    Returns
    -------
    df_targ : pandas.DataFrame
        Dataframe containing merged lithologies/targets
    
    """    
    
    #by default, use the boolean input 
    if target_class=='bool':
        targinterps_df[target_col] = targinterps_df[target_col].where(targinterps_df[target_col]=='1', other='0').astype(int)
        targinterps_df[target_col].fillna(value=0, inplace=True)
    else:
        targinterps_df[target_col].replace('DoNotUse', value=-1, inplace=True)
        targinterps_df[target_col].fillna(value=-2, inplace=True)
        targinterps_df[target_col].astype(np.int8)

    df_targ = pd.merge(df, targinterps_df.set_index(interp_col), right_on=interp_col, left_on='LITHOLOGY', how='left')
    
    return df_targ

#Function to get unique wells

def get_unique_wells(df, wellid_col='API_NUMBER', verbose=False, log=False):
    """Gets unique wells as a dataframe based on a given column name.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing all wells and/or well intervals of interest
    wellid_col : str, default="API_NUMBER"
        Name of column in df containing a unique identifier for each well, by default 'API_NUMBER'. .unique() will be run on this column to get the unique values.
    log : bool, default = False
        Whether to log results to log file

    Returns
    -------
    wellsDF
        DataFrame containing only the unique well IDs
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    #Get Unique well APIs
    uniqueWells = df[wellid_col].unique()
    wellsDF = pd.DataFrame(uniqueWells)
    if verbose:
        print('Number of unique wells in downholeData: '+str(wellsDF.shape[0]))
    wellsDF.columns = ['UNIQUE_ID']
    
    return wellsDF

#Quickly sort dataframe
def sort_dataframe(df, sort_cols=['API_NUMBER','TOP'], remove_nans=True):
    """Function to sort dataframe by one or more columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe to be sorted
    sort_cols : str or list of str, default = ['API_NUMBER','TOP']
        Name(s) of columns by which to sort dataframe, by default ['API_NUMBER','TOP']
    remove_nans : bool, default = True
        Whether or not to remove nans in the process, by default True

    Returns
    -------
    df_sorted : pandas.DataFrame
        Sorted dataframe
    """
    #Sort columns for better processing later
    df_sorted = df.sort_values(sort_cols)
    df_sorted.reset_index(inplace=True, drop=True)
    if remove_nans:
        df_sorted = df_sorted[pd.notna(df_sorted["LITHOLOGY"])]
    return df_sorted