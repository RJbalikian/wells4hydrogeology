import pandas as pd
import numpy as np
import datetime

#The following flags are used to mark the classification method:
#- 0: Not classified
#- 1: Specific Search Term Match
#- 2: wPermits bedrock top pick
#- 3: Intervals >550' below ground surface
#- 4: Wildcard match (startTerm) - no context
#- 5: Wildcard match (startTerm) - more liberal
#- Top of well?

#Define records with full search term
def specificDefine(df, specterms, printouts=False):
    df['FORMATION'] = df['FORMATION'].astype(str)
    specterms['FORMATION'] = specterms['FORMATION'].astype(str)

    print(df.FORMATION)
    df.FORMATION = df.FORMATION.str.casefold()
    specterms['FORMATION'] = specterms['FORMATION'].str.casefold()
    #df['FORMATION'] = df['FORMATION'].str.strip(['.,:?\t\s'])
    #specterms['FORMATION'] = specterms['FORMATION'].str.strip(['.,:?\t\s'])

    print(df.FORMATION)
    specterms.drop_duplicates(subset='FORMATION',keep='last', inplace=True)
    specterms.reset_index(drop=True, inplace=True)
    
    df_Interps = pd.merge(df, specterms.set_index('FORMATION'), on='FORMATION',how='left')
    df_Interps['BEDROCK_FLAG'] = df_Interps['INTERPRETATION'] == 'BEDROCK'
    
    print(df_Interps.FORMATION)
    if printouts:
        print("Records Classified with full search term: "+str(int(df_Interps['CLASS_FLAG'].sum())))
        print("Records Classified with full search term: "+str(round((df_Interps['CLASS_FLAG'].sum()/df_Interps.shape[0])*100,2))+"% of data")
    
    return df_Interps

def splitDefined(df, printouts=False):
    classifedDF= df[df['CLASS_FLAG'].notna()] #Already-classifed data
    searchDF = df[df['CLASS_FLAG'].isna()] #Unclassified data
    
    if printouts:
        print(str(searchDF.shape[0])+' records unclassified; isolated into searchDF.')
    
    return classifedDF, searchDF

#Classify downhole data by the initial substring
def startDefine(df, starterms, printouts=False):
    if printouts:
        estTime = df.shape[0]/3054409 * 6 #It took about 6 minutes to classify data with entire dataframe. This estimates the fraction of that it will take
        nowTime = datetime.datetime.now()
        endTime = nowTime+datetime.timedelta(minutes=estTime)
        print("Start Term process should be done by {:d}:{:02d}".format(endTime.hour, endTime.minute))

    for i,s in enumerate(starterms['FORMATION']):
        df['CLASS_FLAG'].where(~df['FORMATION'].str.startswith(s,na=False),4,inplace=True)
        df['INTERPRETATION'].where(~df['FORMATION'].str.startswith(s,na=False),starterms.loc[i,'INTERPRETATION'],inplace=True)
    df['BEDROCK_FLAG'].loc[df["INTERPRETATION"] == 'BEDROCK']
    
    if printouts:
        print("Records classified with start search term: "+str(int(df['CLASS_FLAG'].count())))
        print("Records classified with start search term: "+str(round((df['CLASS_FLAG'].count()/df.shape[0])*100,2))+"% of remaining data")
        #print("Records classified with both search terms: "+str(round(((df['CLASS_FLAG'].count()+specDF['CLASS_FLAG'].count())/downholeData_Interps.shape[0])*100,2))+"% of all data")
        #This step usually takes about 5-6 minutes
    return df

#Merge data back together
def remergeData(classifieddf, searchdf):
    remergeDF = pd.concat([classifieddf,searchdf], join='inner').sort_index()
    return remergeDF

def depthDefine(dfIN, thresh=550, printouts=False):
    df = dfIN.copy()
    df['CLASS_FLAG'].mask(df['TOP']>thresh, 3 ,inplace=True) #Add a Classification Flag of 3 (bedrock b/c it's deepter than 550') to all records where the top of the interval is >550'
    df['BEDROCK_FLAG'].mask(df['TOP']>thresh, True, inplace=True)

    if printouts:
        if df.CLASS_FLAG.notnull().sum() == 0:
            brDepthClass = 0
        else:
            brDepthClass = df['CLASS_FLAG'].value_counts()[3.0]
        total = dfIN.shape[0]
        print("Records classified as bedrock that were deeper than "+str(thresh)+ "': " + str(brDepthClass))
        print("This represents "+str(round((brDepthClass)*100/total,2))+"% of the unclassified data in this dataframe.")
        
    return df

#Output data that still needs to be defined
def export_toBeDefined(df, outdir):
    #Get directory path correct
    if outdir[-1] != '/' or outdir[-1] != '\\':
        outdir = outdir+'/'
    outdir.replace('\\','/')
    
    todayDate = datetime.date.today()
    todayDateStr = str(todayDate)
    searchDF = df[df['CLASS_FLAG'].isna()]
    
    stillNeededDF=searchDF['FORMATION'].value_counts()
    stillNeededDF.to_csv(outdir+'Stillneed2BeDefined_'+todayDateStr+'.csv')
    return stillNeededDF

def fillUnclassified(df):
    df['CLASS_FLAG'].fillna(0, inplace=True)
    return df

def mergeLithologies(downholedata, targinterps, targetClass='bool'):
    #downholeData = downholedata.copy()
    
    #by default, use the boolean input 
    if targetClass=='bool':
        targinterps['TARGET'] = targinterps['TARGET'].where(targinterps['TARGET']=='1', other='0').astype(int)
        targinterps['TARGET'].fillna(value=0, inplace=True)

    else:
        targinterps['TARGET'].replace('DoNotUse', value=-1, inplace=True)
        targinterps['TARGET'].fillna(value=-2, inplace=True)
        targinterps['TARGET'].astype(np.int8)

    downholeData_targ = pd.merge(downholedata, targinterps.set_index('INTERPRETATION'), right_on='INTERPRETATION',left_on='INTERPRETATION', how='left')
    
    return downholeData_targ

def getUniqueWells(df, wellidCol='API_NUMBER'):
    #Get Unique well APIs
    uniqueWells = df[wellidCol].unique()
    wellsDF = pd.DataFrame(uniqueWells)
    print('Number of unique wells in downholeData: '+str(wellsDF.shape[0]))
    wellsDF.columns = ['UNIQUE_API']
    
    return wellsDF