import pandas as pd
import datetime
import pathlib
import xarray as xr

#Export data 
def export_dataframe(df, procdir, filename, date_stamp=True):
    """Function to export dataframes

    Parameters
    ----------
    df : pandas dataframe, or list of pandas dataframes
        Data frame or list of dataframes to be exported
    procdir : string or pathlib.Path object
        Directory to which to export dataframe object(s) as .csv
    filename : str or list of strings
        Filename(s) of output files
    date_stamp : bool, default=True
        Whether to include a datestamp in the filename. If true, file ends with _yyyy-mm-dd.csv of current date, by default True.
    """

    if date_stamp:
        todayDate = datetime.date.now()
        todayDateStr = '_'+str(todayDate)
    else:
        todayDateStr=''

    if type(procdir) is str or isinstance(procdir, pathlib.PurePath):
        procdir = str(procdir)
        procdir = procdir.replace('\\', '/').replace('\\'[-1], '/')
        if procdir[-1] != '/':
            procdir = procdir + '/'
    else:
        print('Please input string or pathlib object for procdir parameters')
        return

    if type(filename) is str:
        dfOutFile =  procdir+filename+todayDateStr+'.csv'
        df.to_csv(dfOutFile, index_label='ID')
        print('Exported '+filename+todayDateStr+'.csv')
    elif type(filename) is list and type(df) is list and len(df) == len(filename):
        for i, f in enumerate(df):
            fname = filename[i]
            dfOutFile =  procdir+fname+todayDateStr+'.csv'
            f.to_csv(dfOutFile, index_label='ID')
            print('Exported '+fname+todayDateStr+'.csv')

#Export (rio)xarray dataarrays and datasets
def export_grids(grid_data, out_path, filetype='tif', variable_sep=False, date_stamp=True):
    """Function to export grids to files.

    Parameters
    ----------
    grid_data : xarray DataArray or xarray Dataset
        Dataset or dataarray to be exported
    out_path : str or pathlib.Path object
        Output location for data export. If variable_sep=True, this should be a directory. Otherwise, this should also include the filename. The file extension should not be included here.
    filetype : str, optional
        Output filetype. Can either be pickle or any file extension supported by rioxarray.rio.to_raster(). Can either include period or not., by default 'tif'
    variable_sep : bool, optional
        If grid_data is an xarray Dataset, this will export each variable in the dataset as a separate file, including the variable name in the filename, by default False
    date_stamp : bool, optional
        Whether to include a date stamp in the file name., by default True
    """

    #Initialize lists to determine which filetype will be used for export
    ncdfList = ['netcdf', 'ncdf', 'n']
    tifList = ['tif', 'tiff', 'geotiff', 'geotif', 't']
    pickleList = ['pickle', 'pkl', 'p']

    #Format output string(s)
    #Format output filepath
    if type(out_path) is str or isinstance(out_path, pathlib.PurePath):
        if isinstance(out_path, pathlib.PurePath):
            if out_path.parent.exists()==False:
                print('Directory does not exist. Please enter a different value for the out_path parameter.')
                return
        out_path = str(out_path)
        out_path = out_path.replace('\\', '/').replace('\\'[-1], '/')
        if out_path[-1] != '/':
            out_path = out_path + '/'
    else:
        print('Please input string or pathlib object for out_path parameters')
        return
    
    #Format datestamp, if desired in output filename
    if date_stamp:
        todayDate = datetime.date.now()
        todayDateStr = '_'+str(todayDate)
    else:
        todayDateStr=''

    #Ensure the file suffix includes .
    if filetype[0] == '.':
        pass
    else:
        filetype = '.' + filetype

    outPath = out_path+todayDateStr+filetype

    #Do export
    if filetype.lower() in pickleList:
        import pickle
        try:
            with open(outPath, 'wb') as f:
                pickle.dump(grid_data, f)
        except:
            print('An error occured during export.')
            print(outPath, 'could not be exported as a pickle object.')
            print('Try again using different parameters.')
    else:
        import rioxarray as rxr
        try:
            if variable_sep and isinstance(grid_data, xr.Dataset):
                for var in grid_data.data_vars:
                    outPath = out_path+var+todayDateStr+filetype
                    grid_data[var].rio.to_raster(outPath)
            else:
                grid_data.rio.to_raster(outPath)
        except:
            print('An error occured during export.')
            print('{} could not be exported as a {} file.'.format(outPath, filetype))
            print('Try again using different parameters.')

    return