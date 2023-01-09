import pandas as pd
import numpy as np
import pathlib
import datetime
import os

def findMostRecentFiles(dir='../res', globPattern='*'):
    todayDate = datetime.date.today()
    todayDateStr = str(todayDate)

    files = pathlib.Path(dir).rglob(globPattern) #Get all the files that fit the pattern

    fileDates = []
    for f in files: #Get the file dates from their file modification times
        fileDates.append(np.datetime64(datetime.datetime.fromtimestamp(os.path.getmtime(f))))
    globInd = np.argmin(np.datetime64(todayDateStr)- np.array(fileDates)) #Find the index of the most recent file

    #Iterate through glob/files again (need to recreate glob)
    files = pathlib.Path(dir).rglob(globPattern)
    for j, f in enumerate(files):
        if j == globInd:
            mostRecentFile=f
            break
    print('Most Recent version of this file is : '+mostRecentFile.name)
    return mostRecentFile

def filesSetup(db_dir='../res', proc_dir='../out'):
    #Define  filepath variables to be used later for reading/writing files


    rawDirStr = '\\\\isgs-sinkhole\\geophysics\\Balikian\\ISWS_HydroGeo\\RawWellData_OracleDatabase\\TxtData\\'
    processDirStr = '\\\\isgs-sinkhole\\geophysics\\Balikian\\ISWS_HydroGeo\\WellDataAutoClassification\\ProcessedData\\'

    rawDirectory = pathlib.Path(rawDirStr)
    processDirectory = pathlib.Path(processDirStr)

    downholeDataFILE = "lcl_ISGS_DOWNHOLE_DATA.txt"
    headerDataFILE = "lcl_ISGS_HEADER.txt"

    encodeType="latin-1"

    #Define datatypes, to use later
    downholeDataDTYPES = {'ID':np.uint32, "API_NUMBER":np.uint64,"TABLE_NAME":str,"WHO":str,"INTERPRET_DATE":str,"FORMATION":str,"THICKNESS":np.float64,"TOP":np.float64,"BOTTOM":np.float64}
    headerDataDTYPES = {'ID':np.uint32,'API_NUMBER':np.uint64,"TDFORMATION":str,"PRODFORM":str,"TOTAL_DEPTH":np.float64,"SECTION":np.float64,"TWP":np.float64,"TDIR":str,"RNG":np.float64,"RDIR":str,"MERIDIAN":np.float64,"FARM_NAME":str,"NSFOOT":np.float64,"NSDIR":str,"EWFOOT":np.float64,"EWDIR":str,"QUARTERS":str,"ELEVATION":np.float64,"ELEVREF":str,"COMP_DATE":str,"STATUS":str,"FARM_NUM":str,"COUNTY_CODE":np.float64,"PERMIT_NUMBER":str,"COMPANY_NAME":str,"COMPANY_CODE":str,"PERMIT_DATE":str,"CORNER":str,"LATITUDE":np.float64,"LONGITUDE":np.float64,"ENTERED_BY":str,"UPDDATE":str,"ELEVSOURCE":str, "ELEV_FT":np.float64}
