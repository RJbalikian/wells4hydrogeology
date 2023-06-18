import datetime
import logging
import os
import pathlib

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

def run(well_data_file, well_data_columns, 
                     well_metadata_file, well_metadata_columns, 
                     xyz_file):
    """Stil in progress"""
    return