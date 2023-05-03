import datetime
import logging

log_filename = None

#Log data to file
def logger(func):
    def wrapper(*args, **kwargs):
        global log_filename
        log_file = kwargs.pop('log', None)
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
        return result
    return wrapper