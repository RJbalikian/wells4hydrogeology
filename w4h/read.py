import pandas as pd
import numpy as np
import pathlib
import datetime
import os
import json
repoDir = pathlib.Path(os.getcwd())

# Gets the current date for use with in code
def getCurrentDate():
    """ Gets the current date to help with finding the most recent file
        ---------------------
        Parameters:
            None

        ---------------------
        Returns:
            todayDate   : datetime object with today's date
            dateSuffix  : str to use for naming output files
    """
    todayDate = datetime.date.today()
    todayDateStr = str(todayDate)
    dateSuffix = '_'+todayDateStr
    return todayDate, dateSuffix

def findMostRecentFiles(dir=str(repoDir)+'/resources', globPattern='*'):
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

def filesSetup(db_dir=str(repoDir)+'/resources', proc_dir=str(repoDir)+'/out'):
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

#Get filepath of resource in resource folder
def get_resource_path(res):
    repoPath = pathlib.Path(__file__).parent.parent
    repoPathStr = str(repoPath).replace('\\', '/').replace('\\'[0], '/')
    resource = repoPathStr+'/resources/'+res
    return resource

#Read dictionary file into dictionary variable
def read_dict(file='', keytype='np'):
    with open(file, 'r') as f:
        data= f.read()

    jsDict = json.loads(data)
    if keytype=='np':
        for k in jsDict.keys():
            jsDict[k] = getattr(np, jsDict[k])
    
    return jsDict

#Define the datatypes for a dataframe
def defineDataTypes(dfIN, dtypes='', dtypeDir=str(repoDir)+'/resources/', dtypeFile=''):
    df = dfIN.copy()
    
    if dtypeFile != '':
        dtypeFilePATH = pathlib.Path(dtypeDir+dtypeFile)
        dtypes = read_dict(file=dtypeFilePATH)   
    
    for i in range(0, np.shape(df)[1]):
        df.iloc[:,i] = dfIN.iloc[:,i].astype(dtypes[dfIN.iloc[:,i].name])
    return df

def searchTermFilePaths(dictdir=str(repoDir)+'/resources/', specStartPattern='*SearchTerms-Specific*', startGlobPattern = '*SearchTerms-Start*'):
    #Read in dictionary files for downhole data

    #specTermsFile = "SearchTerms-Specific_BedrockOrNo_2022-09.csv" #Specific matches
    #startTermsFile = "SearchTerms-Start_BedrockOrNo.csv" #Wildcard matches for the start of the description
    
    specTermsPATH = findMostRecentFiles(dictdir, specStartPattern)
    startTermsPATH = findMostRecentFiles(dictdir, startGlobPattern)
    
    return specTermsPATH, startTermsPATH

def read_dictionary_terms(dict_file, cols=None, col_types=None, dictionary_type=None, class_flag=1, rem_extra_cols=True):
    #Read files into pandas dataframes
    dict_terms = []
    if type(dict_file) is list:
        for f in dict_file:
            dict_terms.append(pd.read_csv(f))
            if 'ID' in dict_terms.columns:
                dict_terms.set_index('ID', drop=True, inplace=True)
    else:
        dict_terms.append(pd.read_csv(dict_file))
        if 'ID' in dict_terms[-1].columns:
            dict_terms[-1].set_index('ID', drop=True, inplace=True)
        dict_file = [dict_file]

    #Rename important columns
    searchTermList = ['searchterm', 'search', 'exact']
    startTermList = ['startterm', 'start', 'startswith']
                
    #Recast all columns to datatypes of headerData to defined types
    dict_termDtypes = {'FORMATION':str,'INTERPRETATION':str, 'CLASS_FLAG':np.uint8}

    if dictionary_type is None:
        dictionary_type=''

    for i, d in enumerate(dict_terms):
        if dictionary_type.lower() in searchTermList or (dictionary_type=='' and 'spec' in str(dict_file[i]).lower()):
            d['CLASS_FLAG'] = 1
        elif dictionary_type.lower() in startTermList or (dictionary_type=='' and 'start' in str(dict_file[i]).lower()):
            d['CLASS_FLAG'] = 4 #Start term classification flag
        else:
            d['CLASS_FLAG'] = class_flag #Custom classification flag, defined as argument
            #1=specific match,2 (not defined), 3: bedrock classification for obvious bedrock, 4: Wildcard/start term

        if cols is None:
            if 'SearchTerm' in d.columns:
                d.rename(columns={'SearchTerm':'FORMATION'}, inplace=True)
            if 'StartTerm' in d.columns:
                d.rename(columns={'StartTerm':'FORMATION'}, inplace=True)        
            if 'InterpUpdate' in d.columns:
                d.rename(columns={'InterpUpdate':'INTERPRETATION'}, inplace=True)
        else:
            for k in list(cols.keys()):
                d.rename(columns={k:cols[k]}, inplace=True)
        if col_types is None:
            pass
        else:
            for k in list(col_types.keys()):
                d.rename(columns={k:col_types[k]}, inplace=True)
    
        for i in range(0, np.shape(d)[1]):
            if d.iloc[:,i].name in list(dict_termDtypes.keys()):
                d.iloc[:,i] = d.iloc[:,i].astype(dict_termDtypes[d.iloc[:,i].name])
    
        #Delete duplicate definitions
        d.drop_duplicates(subset='FORMATION',inplace=True) #Apparently, there are some duplicate definitions, which need to be deleted first
        d.reset_index(inplace=True, drop=True)

    if len(dict_terms)==1:
        dict_terms = dict_terms[0]
    
    if rem_extra_cols:
        dict_terms = dict_terms[['FORMATION', 'INTERPRETATION', 'CLASS_FLAG']]

    return dict_terms

def readLithologies(lithoDir='', lithFile=''):
    #dictDir = "\\\\isgs-sinkhole\\geophysics\\Balikian\\ISWS_HydroGeo\\WellDataAutoClassification\\SupportingDocs\\"
    if lithoDir =='':
        lithoDir=str(repoDir)+'/resources/'

    if lithFile=='':
        lithFile='Lithology_Interp_FineCoarse.csv'
    
    lithFPath = pathlib.Path(lithoDir+lithFile)
    lithoDF = pd.read_csv(lithFPath, usecols=['LITHOLOGY', 'CODE'])
    lithoDF.rename(columns={'LITHOLOGY':'INTERPRETATION', 'CODE':'TARGET'}, inplace=True)

    return lithoDF