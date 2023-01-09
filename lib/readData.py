import pandas as pd
import numpy as np
import pathlib
import datetime
import os

def getCurrentDate():
    todayDate = datetime.date.today()
    todayDateStr = str(todayDate)
    dateSuffix = '_'+todayDateStr
    return todayDate, dateSuffix

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

    #rawDirStr = '\\\\isgs-sinkhole\\geophysics\\Balikian\\ISWS_HydroGeo\\RawWellData_OracleDatabase\\TxtData\\'
    #processDirStr = '\\\\isgs-sinkhole\\geophysics\\Balikian\\ISWS_HydroGeo\\WellDataAutoClassification\\ProcessedData\\'

    rawDirectory = pathlib.Path(db_dir)
    processDirectory = pathlib.Path(proc_dir)

    downholeDataFILE = findMostRecentFiles(db_dir, '*ISGS_DOWNHOLE_DATA*.txt')
    headerDataFILE = findMostRecentFiles(db_dir, '*ISGS_HEADER*.txt')
    xyzInFILE = findMostRecentFiles(db_dir, '*xyzData*')

    downholeDataPATH = pathlib.path(downholeDataFILE)
    headerDataPATH = pathlib.path(headerDataFILE)
    xyzInPATH = pathlib.path(xyzInFILE)

    encodeType="latin-1"

    print('Using the following files:\n')
    print(downholeDataFILE)
    print(headerDataFILE)
    print(xyzInFILE)

    #Define datatypes, to use later
    #downholeDataDTYPES = {'ID':np.uint32, "API_NUMBER":np.uint64,"TABLE_NAME":str,"WHO":str,"INTERPRET_DATE":str,"FORMATION":str,"THICKNESS":np.float64,"TOP":np.float64,"BOTTOM":np.float64}
    #headerDataDTYPES = {'ID':np.uint32,'API_NUMBER':np.uint64,"TDFORMATION":str,"PRODFORM":str,"TOTAL_DEPTH":np.float64,"SECTION":np.float64,"TWP":np.float64,"TDIR":str,"RNG":np.float64,"RDIR":str,"MERIDIAN":np.float64,"FARM_NAME":str,"NSFOOT":np.float64,"NSDIR":str,"EWFOOT":np.float64,"EWDIR":str,"QUARTERS":str,"ELEVATION":np.float64,"ELEVREF":str,"COMP_DATE":str,"STATUS":str,"FARM_NUM":str,"COUNTY_CODE":np.float64,"PERMIT_NUMBER":str,"COMPANY_NAME":str,"COMPANY_CODE":str,"PERMIT_DATE":str,"CORNER":str,"LATITUDE":np.float64,"LONGITUDE":np.float64,"ENTERED_BY":str,"UPDDATE":str,"ELEVSOURCE":str, "ELEV_FT":np.float64}
    return downholeDataPATH, headerDataPATH, xyzInPATH

def readRawTxtData(rawdir='../res', downholefile='', headerfile='', encoding='latin-1'):
    
    headers_useCols = ['API_NUMBER',"TOTAL_DEPTH","SECTION","TWP","TDIR","RNG","RDIR","MERIDIAN","QUARTERS","ELEVATION","ELEVREF","COUNTY_CODE","LATITUDE","LONGITUDE","ELEVSOURCE"]
    downhole_useCols = ["API_NUMBER","TABLE_NAME","FORMATION","THICKNESS","TOP","BOTTOM"]
    
    encodeType=encoding
    
    downholeDataIN = pd.read_csv(rawdir+downholefile, sep=',', header='infer', encoding=encodeType, usecols=downhole_useCols)
    headerDataIN = pd.read_csv(rawdir+headerfile, sep=',', header='infer', encoding=encodeType, usecols=headers_useCols)
        
    downholeDataIN = downholeDataIN.dropna(subset=['API_NUMBER']) #Drop data with no API
    
    #Drop data with no or missing location information
    headerDataIN = headerDataIN.dropna(subset=['LATITUDE']) 
    headerDataIN = headerDataIN.dropna(subset=['LONGITUDE'])
    
    #Reset index so index goes from 0 in numerical/integer order
    headerDataIN.reset_index(inplace=True, drop=True)
    downholeDataIN.reset_index(inplace=True, drop=True)
    
    print('Downhole Data has ' + str(downholeDataIN.shape[0])+' valid well records.')
    print("Header Data has "+str(headerDataIN.shape[0])+" unique wells with valid location information.")
    
    return headerDataIN, downholeDataIN

#Read file with xyz data
def readXYZData(rawdir='../res', xyzfile=''):
    xyzDTypes = {'ID':np.uint32,'API_NUMBER':np.uint64,'LATITUDE':np.float64,'LONGITUDE':np.float64,'ELEV_FT':np.float64}
    xyzDataIN = pd.read_csv(rawdir+xyzfile, sep=',', header='infer', dtype=xyzDTypes, index_col='ID')
    
    return xyzDataIN

def defineDataTypes(dfIN, dtypes):
    df = dfIN.copy()
    
    for i in range(0, np.shape(df)[1]):
        df.iloc[:,i] = dfIN.iloc[:,i].astype(dtypes[dfIN.iloc[:,i].name])
    return df