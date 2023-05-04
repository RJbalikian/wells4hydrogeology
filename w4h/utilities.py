import datetime
import logging

log_filename = None

#Log data to file
def logger(func):
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
    return wrapper