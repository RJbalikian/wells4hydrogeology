import pandas as pd
import numpy as np
import pathlib
import datetime
import os
import json
repoDir = pathlib.Path(os.getcwd())

def getCurrentDate():
    todayDate = datetime.date.today()
    todayDateStr = str(todayDate)
    dateSuffix = '_'+todayDateStr
    return todayDate, dateSuffix

def findMostRecentFiles(dir=str(repoDir)+'/res', globPattern='*'):
    todayDate = datetime.date.today()
    todayDateStr = str(todayDate)

    files = pathlib.Path(dir).rglob(globPattern) #Get all the files that fit the pattern

    fileDates = []
    for f in files: #Get the file dates from their file modification times
        fileDates.append(np.datetime64(datetime.datetime.fromtimestamp(os.path.getmtime(f))))
    globInd = np.argmin(np.datetime64(todayDateStr) - np.array(fileDates)) #Find the index of the most recent file

    #Iterate through glob/files again (need to recreate glob)
    files = pathlib.Path(dir).rglob(globPattern)
    for j, f in enumerate(files):
        if j == globInd:
            mostRecentFile=f
            break
    print('Most Recent version of this file is : '+mostRecentFile.name)
    return mostRecentFile

def filesSetup(db_dir=str(repoDir)+'/res', proc_dir=str(repoDir)+'/out'):
    #Define  filepath variables to be used later for reading/writing files

    #rawDirStr = '\\\\isgs-sinkhole\\geophysics\\Balikian\\ISWS_HydroGeo\\RawWellData_OracleDatabase\\TxtData\\'
    #processDirStr = '\\\\isgs-sinkhole\\geophysics\\Balikian\\ISWS_HydroGeo\\WellDataAutoClassification\\ProcessedData\\'

    rawDirectory = pathlib.Path(db_dir)
    processDirectory = pathlib.Path(proc_dir)

    downholeDataFILE = findMostRecentFiles(db_dir, '*ISGS_DOWNHOLE_DATA*.txt')
    headerDataFILE = findMostRecentFiles(db_dir, '*ISGS_HEADER*.txt')
    xyzInFILE = findMostRecentFiles(db_dir, '*xyzData*')

    downholeDataPATH = pathlib.Path(downholeDataFILE)
    headerDataPATH = pathlib.Path(headerDataFILE)
    xyzInPATH = pathlib.Path(xyzInFILE)

    encodeType="latin-1"

    print('Using the following files:\n')
    print(downholeDataFILE)
    print(headerDataFILE)
    print(xyzInFILE)

    #Define datatypes, to use later
    #downholeDataDTYPES = {'ID':np.uint32, "API_NUMBER":np.uint64,"TABLE_NAME":str,"WHO":str,"INTERPRET_DATE":str,"FORMATION":str,"THICKNESS":np.float64,"TOP":np.float64,"BOTTOM":np.float64}
    #headerDataDTYPES = {'ID':np.uint32,'API_NUMBER':np.uint64,"TDFORMATION":str,"PRODFORM":str,"TOTAL_DEPTH":np.float64,"SECTION":np.float64,"TWP":np.float64,"TDIR":str,"RNG":np.float64,"RDIR":str,"MERIDIAN":np.float64,"FARM_NAME":str,"NSFOOT":np.float64,"NSDIR":str,"EWFOOT":np.float64,"EWDIR":str,"QUARTERS":str,"ELEVATION":np.float64,"ELEVREF":str,"COMP_DATE":str,"STATUS":str,"FARM_NUM":str,"COUNTY_CODE":np.float64,"PERMIT_NUMBER":str,"COMPANY_NAME":str,"COMPANY_CODE":str,"PERMIT_DATE":str,"CORNER":str,"LATITUDE":np.float64,"LONGITUDE":np.float64,"ENTERED_BY":str,"UPDDATE":str,"ELEVSOURCE":str, "ELEV_FT":np.float64}
    return downholeDataPATH, headerDataPATH, xyzInPATH

def readRawTxtData(rawdir='', downholefile='', headerfile='', encoding='latin-1'):
    
    headers_useCols = ['API_NUMBER',"TOTAL_DEPTH","SECTION","TWP","TDIR","RNG","RDIR","MERIDIAN","QUARTERS","ELEVATION","ELEVREF","COUNTY_CODE","LATITUDE","LONGITUDE","ELEVSOURCE"]
    downhole_useCols = ["API_NUMBER","TABLE_NAME","FORMATION","THICKNESS","TOP","BOTTOM"]
    
    encodeType=encoding
    
    downholeDataIN = pd.read_csv(rawdir+str(downholefile), sep=',', header='infer', encoding=encodeType, usecols=downhole_useCols)
    headerDataIN = pd.read_csv(rawdir+str(headerfile), sep=',', header='infer', encoding=encodeType, usecols=headers_useCols)
        
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
def readXYZData(rawdir='', xyzfile=''):
    xyzDTypes = {'ID':np.uint32,'API_NUMBER':np.uint64,'LATITUDE':np.float64,'LONGITUDE':np.float64,'ELEV_FT':np.float64}
    xyzDataIN = pd.read_csv(rawdir+str(xyzfile), sep=',', header='infer', dtype=xyzDTypes, index_col='ID')
    
    return xyzDataIN

#Only keep 'essential' columns to decrease size of dataframes/speed up processing
#def essentialCols(downholedata, headerdata):
#    downhole_keepCols = ["API_NUMBER","TABLE_NAME","FORMATION","THICKNESS","TOP","BOTTOM"]
#    headers_keepCols = ['API_NUMBER',"TOTAL_DEPTH","SECTION","TWP","TDIR","RNG","RDIR","MERIDIAN","QUARTERS","ELEVATION","ELEVREF","COUNTY_CODE","LATITUDE","LONGITUDE","ELEVSOURCE"]
    
#    downholeData = downholedata[downhole_keepCols]
#    headerData = headerdata[headers_keepCols]
    
#    return downholeData, headerData

def readDataTypeDict(file=''):
    with open(file, 'r') as f:
        data= f.read()

    jsDict = json.loads(data)
    for k in jsDict.keys():
        jsDict[k] = getattr(np, jsDict[k])
        
    return jsDict

#Define the datatypes for a dataframe
def defineDataTypes(dfIN, dtypes='', dtypeDir=str(repoDir)+'/res/', dtypeFile=''):
    df = dfIN.copy()
    
    if dtypeFile != '':
        dtypeFilePATH = pathlib.Path(dtypeDir+dtypeFile)
        dtypes = readDataTypeDict(file=dtypeFilePATH)   
    
    for i in range(0, np.shape(df)[1]):
        df.iloc[:,i] = dfIN.iloc[:,i].astype(dtypes[dfIN.iloc[:,i].name])
    return df

def searchTermFilePaths(dictdir=str(repoDir)+'/res/', specStartPattern='*SearchTerms-Specific*', startGlobPattern = '*SearchTerms-Start*'):
    #Read in dictionary files for downhole data

    #specTermsFile = "SearchTerms-Specific_BedrockOrNo_2022-09.csv" #Specific matches
    #startTermsFile = "SearchTerms-Start_BedrockOrNo.csv" #Wildcard matches for the start of the description
    
    specTermsPATH = findMostRecentFiles(dictdir, specStartPattern)
    startTermsPATH = findMostRecentFiles(dictdir, startGlobPattern)
    
    return specTermsPATH, startTermsPATH

def readSearchTerms(specfile, startfile, dictdir=str(repoDir)+'/res/'):
    #Read files into pandas dataframes
    specTerms = pd.read_csv(specfile, index_col='ID')
    startTerms = pd.read_csv(startfile, index_col='ID')

    #Rename important columns
    specTerms.rename(columns={'SearchTerm':'FORMATION', 'InterpUpdate':'INTERPRETATION'}, inplace=True)
    specTerms['CLASS_FLAG'] = 1
    startTerms.rename(columns={'StartTerm':'FORMATION', 'InterpUpdate':'INTERPRETATION'}, inplace=True)
    startTerms['CLASS_FLAG'] = 4

    #Recast all columns to datatypes of headerData to defined types
    specTermsDtypes = {'FORMATION':str,'INTERPRETATION':str, 'CLASS_FLAG':np.uint8}
    startTermsDtypes = {'FORMATION':str,'INTERPRETATION':str, 'CLASS_FLAG':np.uint8}

    for i in range(0, np.shape(specTerms)[1]):
        specTerms.iloc[:,i] = specTerms.iloc[:,i].astype(specTermsDtypes[specTerms.iloc[:,i].name])
        startTerms.iloc[:,i] = startTerms.iloc[:,i].astype(startTermsDtypes[startTerms.iloc[:,i].name])
    
    #Delete duplicate definitions
    specTerms.drop_duplicates(subset='FORMATION',inplace=True) #Apparently, there are some duplicate definitions, which need to be deleted first
    startTerms.drop_duplicates(subset='FORMATION',inplace=True) #Apparently, there are some duplicate definitions, which need to be deleted first
    specTerms.reset_index(inplace=True, drop=True)
    startTerms.reset_index(inplace=True, drop=True)
    
    return specTerms, startTerms

def readLithologies(lithoDir='', lithFile=''):
    #dictDir = "\\\\isgs-sinkhole\\geophysics\\Balikian\\ISWS_HydroGeo\\WellDataAutoClassification\\SupportingDocs\\"
    if lithoDir =='':
        lithoDir=str(repoDir)+'/res/'

    if lithFile=='':
        lithFile='Lithology_Interp_FineCoarse.csv'
    
    lithFPath = pathlib.Path(lithoDir+lithFile)
    lithoDF = pd.read_csv(lithFPath, usecols=['LITHOLOGY', 'CODE'])
    lithoDF['CODE'] = lithoDF['CODE'].where(lithoDF['CODE']=='1', other=0).astype(int)
    lithoDF.rename(columns={'LITHOLOGY':'INTERPRETATION', 'CODE':'TARGET'}, inplace=True)

    return lithoDF