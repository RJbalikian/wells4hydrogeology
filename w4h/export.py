import datetime
import inspect
import pathlib

import pandas as pd
import xarray as xr

from w4h import logger_function

#Export data

def export_dataframe(df, out_dir, filename, date_stamp=True, log=False):
    """Function to export dataframes

    Parameters
    ----------
    df : pandas dataframe, or list of pandas dataframes
        Data frame or list of dataframes to be exported
    out_dir : string or pathlib.Path object
        Directory to which to export dataframe object(s) as .csv
    filename : str or list of strings
        Filename(s) of output files
    date_stamp : bool, default=True
        Whether to include a datestamp in the filename. If true, file ends with _yyyy-mm-dd.csv of current date, by default True.
    log : bool, default = True
        Whether to log inputs and outputs to log file.        
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    if date_stamp:
        todayDate = datetime.date.today()
        todayDateStr = '_'+str(todayDate)
    else:
        todayDateStr=''

    if type(out_dir) is str or isinstance(out_dir, pathlib.PurePath):
        out_dir = str(out_dir)
        out_dir = out_dir.replace('\\', '/').replace('\\'[-1], '/')
        if out_dir[-1] != '/':
            out_dir = out_dir + '/'
    else:
        print('Please input string or pathlib object for out_dir parameters')
        return

    if type(filename) is str:
        dfOutFile =  out_dir+filename+todayDateStr+'.csv'
        df.to_csv(dfOutFile, index_label='ID')
        print('Exported '+filename+todayDateStr+'.csv')
    elif type(filename) is list and type(df) is list and len(df) == len(filename):
        for i, f in enumerate(df):
            fname = filename[i]
            dfOutFile =  out_dir+fname+todayDateStr+'.csv'
            f.to_csv(dfOutFile, index_label='ID')
            print('Exported '+fname+todayDateStr+'.csv')

#Export (rio)xarray dataarrays and datasets

def export_grids(grid_data, out_path, file_id='',filetype='tif', variable_sep=True, date_stamp=True, verbose=False, log=False):
    """Function to export grids to files.

    Parameters
    ----------
    grid_data : xarray DataArray or xarray Dataset
        Dataset or dataarray to be exported
    out_path : str or pathlib.Path object
        Output location for data export. If variable_sep=True, this should be a directory. Otherwise, this should also include the filename. The file extension should not be included here.
    file_id : str, optional
        If specified, will add this after 'LayerXX' or 'AllLayers' in the filename, just before datestamp, if used. Example filename for file_id='Coarse': Layer1_Coarse_2023-04-18.tif.
    filetype : str, optional
        Output filetype. Can either be pickle or any file extension supported by rioxarray.rio.to_raster(). Can either include period or not., by default 'tif'
    variable_sep : bool, optional
        If grid_data is an xarray Dataset, this will export each variable in the dataset as a separate file, including the variable name in the filename, by default False
    date_stamp : bool, optional
        Whether to include a date stamp in the file name., by default True
    log : bool, default = True
        Whether to log inputs and outputs to log file.        
    """
    logger_function(log, locals(), inspect.currentframe().f_code.co_name)

    #Initialize lists to determine which filetype will be used for export
    ncdfList = ['netcdf', 'ncdf', 'n']
    tifList = ['tif', 'tiff', 'geotiff', 'geotif', 't']
    pickleList = ['pickle', 'pkl', 'p']

    filenames = []

    #Format output string(s)
    #Format output filepath
    if type(out_path) is str or isinstance(out_path, pathlib.PurePath):
        if isinstance(out_path, pathlib.PurePath):
            pass
        else:
            out_path = pathlib.Path(out_path)
            
        if out_path.parent.exists() == False:
            print('Directory does not exist. Please enter a different value for the out_path parameter.')
            return        

        if out_path.is_dir():
            if isinstance(grid_data, xr.DataArray):
                if variable_sep:
                    lyrs = grid_data.coords['Layer'].values
                    filenames = []
                    for l in lyrs:
                        filenames.append('Layer'+str(l))
                else:
                    filenames = ['AllLayers']
            if isinstance(grid_data, xr.Dataset):
                if variable_sep:
                    filenames = []
                    for var in grid_data:
                        filenames.append(var)
                else:
                    filenames = ['AllLayers']    
        else:
            filenames = [out_path.stem]
            out_path = out_path.parent

    else:
        print('Please input string or pathlib object for out_path parameters')
        return
    
    #Format datestamp, if desired in output filename
    if date_stamp:
        todayDate = datetime.date.today()
        todayDateStr = '_'+str(todayDate)
    else:
        todayDateStr=''

    #Ensure the file suffix includes .
    if filetype[0] == '.':
        pass
    else:
        filetype = '.' + filetype

    if file_id != '':
        file_id = '_'+file_id

    out_path = out_path.as_posix()+'/'
    
    if verbose:
        print('Export filepath(s):')
        
    outPaths = []
    for f in filenames:
        currOutPath = out_path+f+file_id+todayDateStr+filetype
        outPaths.append(currOutPath)
        if verbose:
            print('\t {}'.format(currOutPath))
        
    #Do export
    if filetype.lower() in pickleList:
        import pickle
        for op in outPaths:
            try:
                with open(op, 'wb') as f:
                    pickle.dump(grid_data, f)
            except:
                print('An error occured during export.')
                print(op, 'could not be exported as a pickle object.')
                print('Try again using different parameters.')
    else:
        import rioxarray as rxr
        try:
            if isinstance(grid_data, xr.Dataset):
                if variable_sep:
                    for i, var in enumerate(grid_data.data_vars):
                        grid_data[var].rio.to_raster(outPaths[i])
                else:
                    grid_data.rio.to_raster(outPaths[0])
            elif isinstance(grid_data, xr.DataArray):
                if variable_sep:
                    lyrs = grid_data.coords['Layer'].values
                    for i, l in enumerate(lyrs):
                        out_grid = grid_data.sel(Layer = l).copy()
                        out_grid.rio.to_raster(outPaths[i])
                else:
                    grid_data.rio.to_raster(outPaths[0])
            else:
                grid_data.rio.to_raster(outPaths[0])
        except:
            print('An error occured during export.')
            print('{} could not be exported as {} file.'.format(outPaths, filetype))
            print('Try again using different parameters.')

    return
