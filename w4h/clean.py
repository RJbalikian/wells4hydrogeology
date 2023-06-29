import inspect

import numpy as np
import pandas as pd

from w4h import logger_function

#This function removes all data from the downholeData table where there is no location information (in the headerData table). This includes elevation info too
def remove_nonlocated(df, metadata_df, verbose=False, log=False):
    """Function to remove wells and well intervals where there is no location information

    Parameters
    ----------
    df : pandas.DataFrame
        Pandas dataframe containing well descriptions
    metadata_DF : pandas.DataFrame
        Pandas dataframe containing metadata, including well locations (e.g., Latitude/Longitude)
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    df : pandas.DataFrame
        Pandas dataframe containing only data with location information
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    before = df.shape[0] #Extract length of data before this process

    #Create Merged dataset only with data where wells exist in both databases (i.e., well has data and location info)
    df = pd.merge(df, metadata_df.set_index('API_NUMBER'), on='API_NUMBER', how='left', indicator='Exist')
    df['Existbool'] = np.where(df['Exist'] == 'both', True, False)
    df = df[df['Existbool']==True].drop(['Exist','Existbool'], axis=1)
    
    #Create new downhole data table with only relevant records and columns
    keepCols=['API_NUMBER','TABLE_NAME','FORMATION','THICKNESS','TOP','BOTTOM']
    df = df[keepCols].copy()
    if verbose:
        after = df.shape[0]
        print(str(before-after)+' records removed without location information.')
        print(str(df.shape[0])+' wells remain from '+str(df['API_NUMBER'].unique().shape[0])+' geolocated wells in study area.')
    return df

#Function to remove data (intended for headerData) without surface topography information
##THIS ASSUMES AND SHOULD ONLY BE RUN AFTER ALL DESIRED SURFACE TOPO DATASETS HAVE BEEN MERGED/ADDED
def remove_no_topo(df, zcol='ELEV_FT', no_data_val='', verbose=False, log=False):
    """Function to remove wells that do not have topography data (needed for layer selection later).

    This function is intended to be run on the metadata table after elevations have attempted to been added.

    Parameters
    ----------
    df : pandas.DataFrame
        Pandas dataframe containing elevation information.
    zcol : str
        Name of elevation column
    no_data_val : any
        Value in dataset that indicates no data is present (replaced with np.nan)
    verbose : bool, optional
        Whether to print outputs, by default True
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    pandas.DataFrame
        Pandas dataframe with intervals with no topography removed.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    before = df.shape[0]
    
    df[zcol].replace(no_data_val, np.nan, inplace=True)
    df.dropna(subset=[zcol], inplace=True)
    
    if verbose:
        after = df.shape[0]
        print('Well records removed: '+str(before-after))
        print("Number of rows before dropping those without surface elevation information: "+str(before))
        print("Number of rows after dropping those without surface elevation information: "+str(after))
    
    return df

#This function drops all records in the downholedata with no depth information (either top or bottom depth of well interval)
def remove_no_depth(df, top_col='TOP', bottom_col='BOTTOM', no_data_val='', verbose=False, log=False):
    """Function to remove well intervals with no depth information

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing well descriptions
    top_col : str, optional
        Name of column containing information on the top of the well intervals, by default 'TOP'
    bottom_col : str, optional
        Name of column containing information on the bottom of the well intervals, by default 'BOTTOM'
    no_data_val : any, optional
        No data value in the input data, used by this function to indicate that depth data is not there, to be replaced by np.nan, by default ''
    verbose : bool, optional
        Whether to print results to console, by default False
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    df : pandas.DataFrame
        Dataframe with depths dropped
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    #Replace empty cells in top and bottom columns with nan
    df[top_col] = df[top_col].replace(no_data_val, np.nan)
    df[bottom_col] = df[bottom_col].replace(no_data_val, np.nan)
    
    #Calculate number of rows before dropping
    before = df.shape[0]

    #Drop records without depth information
    df = df.dropna(subset=[top_col])
    df = df.dropna(subset=[bottom_col])
    df.reset_index(inplace=True, drop=True) #Reset index
  
    if verbose:
        print("Number of rows before dropping those without record depth information: " + str(before))
        print("Number of rows after dropping those without record depth information: " + str(df.shape[0]))
        print('Number of well records without formation information deleted: ' + str(before - df.shape[0]))
    
    return df

#This function drops all records in downholeData with bad depth information (where the bottom of a record is nearer to the surface than the top)
def remove_bad_depth(df, top_col='TOP', bottom_col='BOTTOM', depth_type='depth', verbose=False, log=False):
    """Function to remove all records in the dataframe with well interpretations where the depth information is bad (i.e., where the bottom of the record is neerer to the surface than the top)

    Parameters
    ----------
    df : pandas.DataFrame
        Pandas dataframe containing the well records and descriptions for each interval
    top_col : str, default='TOP'
        The name of the column containing the depth or elevation for the top of the interval, by default 'TOP'
    bottom_col : str, default='BOTTOM'
        The name of the column containing the depth or elevation for the bottom of each interval, by default 'BOTTOM'
    depth_type : str, {'depth', 'elevation'}
        Whether the table is organized by depth or elevation. If depth, the top column will have smaller values than the bottom column. If elevation, the top column will have higher values than the bottom column, by default 'depth'
    verbose : bool, default = False
        Whether to print results to the terminal, by default False
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    pandas.Dataframe
        Pandas dataframe with the records remvoed where the top is indicatd to be below the bottom.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if depth_type.lower() =='depth':
        df['THICKNESS'] = df[bottom_col] - df[top_col] #Calculate interval thickness
    elif depth_type.lower() =='elevation' or depth_type=='elev':
        df['THICKNESS'] = df[top_col] - df[bottom_col] #Calculate interval thickness
    before = df.shape[0] #Calculate number of rows before dropping
    df = df[(df['THICKNESS'] >= 0)] #Only include rows where interval thickness is positive (bottom is deeper than top)
    df.reset_index(inplace=True, drop=True) #Reset index

    if verbose:
        print("Number of rows before dropping those with obviously bad depth information: "+str(before))
        print("Number of rows after dropping those with obviously bad depth information: "+str(df.shape[0]))
        print('Well records deleted: '+str(before-df.shape[0]))
    return df

#This function drops all records in downholeData with no formation in formation in the description fiel
def remove_no_formation(df, description_col='FORMATION', no_data_val='', verbose=False, log=False):
    """Function that removes all records in the dataframe containing the well descriptions where no description is given.

    Parameters
    ----------
    df : pandas.DataFrame
        Pandas dataframe containing the well records with their individual descriptions
    description_col : str, optional
        Name of the column containing the geologic description of each interval, by default 'FORMATION'
    no_data_val : str, optional
        The value expected if the column is empty or there is no data. These will be replaced by np.nan before being removed, by default ''
    verbose : bool, optional
        Whether to print the results of this step to the terminal, by default False
    log : bool, default = False
        Whether to log results to log file, by default False
        
    Returns
    -------
    pandas.DataFrame
        Pandas dataframe with records with no description removed.
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    #Replace empty cells in formation column with nans
    df[description_col] = df[description_col].replace(no_data_val, np.nan) 
    before = df.shape[0] #Calculate number of rows before dropping

    #Drop records without FORMATION information
    df = df.dropna(subset=[description_col])
    df.reset_index(inplace=True, drop=True) #Reset index

    if verbose:
        print("Number of rows before dropping those without FORMATION information: "+str(before))
        print("Number of rows after dropping those without FORMATION information: "+str(df.shape[0]))
        print('Well records deleted: '+str(before-df.shape[0]))
        
    return df
