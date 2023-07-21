import inspect

import numpy as np
import pandas as pd

from w4h import logger_function

#This function removes all data from the downholeData table where there is no location information (in the headerData table). This includes elevation info too
def remove_nonlocated(df_with_locations, xcol='LONGITUDE', ycol='LATITUDE', no_data_val_table='', verbose=False, log=False):
    """Function to remove wells and well intervals where there is no location information

    Parameters
    ----------
    df_with_locations : pandas.DataFrame
        Pandas dataframe containing well descriptions
    metadata_DF : pandas.DataFrame
        Pandas dataframe containing metadata, including well locations (e.g., Latitude/Longitude)
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    df_with_locations : pandas.DataFrame
        Pandas dataframe containing only data with location information
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    before = df_with_locations.shape[0] #Extract length of data before this process

    df_with_locations[xcol].replace(no_data_val_table, np.nan, inplace=True)
    df_with_locations[ycol].replace(no_data_val_table, np.nan, inplace=True)
    
    df_with_locations.dropna(subset=xcol, inplace=True)
    df_with_locations.dropna(subset=ycol, inplace=True)
    
    if verbose:
        after = df_with_locations.shape[0]
        print('Removed well records with no location information. ')
        print("\tNumber of records before removing: "+str(before))
        print("\tNumber of records after removing: "+str(after))
        print("\t\t{} wells records removed without location information".format(before-after))

    return df_with_locations

#Function to remove data (intended for headerData) without surface topography information
##THIS ASSUMES AND SHOULD ONLY BE RUN AFTER ALL DESIRED SURFACE TOPO DATASETS HAVE BEEN MERGED/ADDED
def remove_no_topo(df_with_topo, zcol='ELEVATION', no_data_val_table='', verbose=False, log=False):
    """Function to remove wells that do not have topography data (needed for layer selection later).

    This function is intended to be run on the metadata table after elevations have attempted to been added.

    Parameters
    ----------
    df_with_topo : pandas.DataFrame
        Pandas dataframe containing elevation information.
    zcol : str
        Name of elevation column
    no_data_val_table : any
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

    before = df_with_topo.shape[0]
    
    df_with_topo[zcol].replace(no_data_val_table, np.nan, inplace=True)
    df_with_topo.dropna(subset=[zcol], inplace=True)
    
    if verbose:
        after = df_with_topo.shape[0]
        print('Removed well records with no surface elevation information. ')
        print("\tNumber of records before removing: "+str(before))
        print("\tNumber of records after removing: "+str(after))
        print("\t\t{} wells records removed without surface elevation information".format(before-after))
    
    return df_with_topo

#This function drops all records in the downholedata with no depth information (either top or bottom depth of well interval)
def remove_no_depth(df_with_depth, top_col='TOP', bottom_col='BOTTOM', no_data_val_table='', verbose=False, log=False):
    """Function to remove well intervals with no depth information

    Parameters
    ----------
    df_with_depth : pandas.DataFrame
        Dataframe containing well descriptions
    top_col : str, optional
        Name of column containing information on the top of the well intervals, by default 'TOP'
    bottom_col : str, optional
        Name of column containing information on the bottom of the well intervals, by default 'BOTTOM'
    no_data_val_table : any, optional
        No data value in the input data, used by this function to indicate that depth data is not there, to be replaced by np.nan, by default ''
    verbose : bool, optional
        Whether to print results to console, by default False
    log : bool, default = False
        Whether to log results to log file, by default False

    Returns
    -------
    df_with_depth : pandas.DataFrame
        Dataframe with depths dropped
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    #Replace empty cells in top and bottom columns with nan
    df_with_depth[top_col] = df_with_depth[top_col].replace(no_data_val_table, np.nan)
    df_with_depth[bottom_col] = df_with_depth[bottom_col].replace(no_data_val_table, np.nan)
    
    #Calculate number of rows before dropping
    before = df_with_depth.shape[0]

    #Drop records without depth information
    df_with_depth = df_with_depth.dropna(subset=[top_col])
    df_with_depth = df_with_depth.dropna(subset=[bottom_col])
    df_with_depth.reset_index(inplace=True, drop=True) #Reset index
  
    if verbose:
        after = df_with_depth.shape[0]
        print('Removed well records with no depth information. ')
        print("\tNumber of records before removing: "+str(before))
        print("\tNumber of records after removing: "+str(after))
        print("\t\t{} well records removed without depth information".format(before-after))
    
    return df_with_depth

#This function drops all records in downholeData with bad depth information (where the bottom of a record is nearer to the surface than the top)
def remove_bad_depth(df_with_depth, top_col='TOP', bottom_col='BOTTOM', depth_type='depth', verbose=False, log=False):
    """Function to remove all records in the dataframe with well interpretations where the depth information is bad (i.e., where the bottom of the record is neerer to the surface than the top)

    Parameters
    ----------
    df_with_depth : pandas.DataFrame
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
        df_with_depth['THICKNESS'] = df_with_depth[bottom_col] - df_with_depth[top_col] #Calculate interval thickness
    elif depth_type.lower() =='elevation' or depth_type=='elev':
        df_with_depth['THICKNESS'] = df_with_depth[top_col] - df_with_depth[bottom_col] #Calculate interval thickness
    before = df_with_depth.shape[0] #Calculate number of rows before dropping
    df_with_depth = df_with_depth[(df_with_depth['THICKNESS'] >= 0)] #Only include rows where interval thickness is positive (bottom is deeper than top)
    df_with_depth.reset_index(inplace=True, drop=True) #Reset index

    if verbose:
        after = df_with_depth.shape[0]
        print('Removed well records with obviously bad depth information. ')
        print("\tNumber of records before removing: "+str(before))
        print("\tNumber of records after removing: "+str(after))
        print("\t\t{} well records removed without depth information".format(before-after))

    return df_with_depth

#This function drops all records in downholeData with no formation in formation in the description fiel
def remove_no_description(df_with_descriptions, description_col='FORMATION', no_data_val_table='', verbose=False, log=False):
    """Function that removes all records in the dataframe containing the well descriptions where no description is given.

    Parameters
    ----------
    df_with_descriptions : pandas.DataFrame
        Pandas dataframe containing the well records with their individual descriptions
    description_col : str, optional
        Name of the column containing the geologic description of each interval, by default 'FORMATION'
    no_data_val_table : str, optional
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
    df_with_descriptions[description_col] = df_with_descriptions[description_col].replace(no_data_val_table, np.nan) 
    before = df_with_descriptions.shape[0] #Calculate number of rows before dropping

    #Drop records without FORMATION information
    df_with_descriptions = df_with_descriptions.dropna(subset=[description_col])
    df_with_descriptions.reset_index(inplace=True, drop=True) #Reset index

    if verbose:
        after = df_with_descriptions.shape[0]
        print('Removed well records without geologic descriptions. ')
        print("\tNumber of records before removing: "+str(before))
        print("\tNumber of records after removing: "+str(after))
        print("\t\t{} well records removed without geologic descriptions".format(before-after))

    return df_with_descriptions
