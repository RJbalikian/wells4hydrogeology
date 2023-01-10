import numpy as np
import pandas as pd

#This function removes all data from the downholeData table where there is no location information (in the headerData table). This includes elevation info too
def removeNonlocatedData(downholeData, headerData):
    before = downholeData.shape[0] #Extract length of data before this process

    #Create Merged dataset only with data where wells exist in both databases (i.e., well has data and location info)
    df = pd.merge(downholeData, headerData.set_index('API_NUMBER'), on='API_NUMBER', how='left', indicator='Exist')
    df['Existbool'] = np.where(df['Exist'] == 'both', True, False)
    df = df[df['Existbool']==True].drop(['Exist','Existbool'], axis=1)
    
    #Create new downhole data table with only relevant records and columns
    keepCols=['API_NUMBER','TABLE_NAME','FORMATION','THICKNESS','TOP','BOTTOM']
    downholeData = df[keepCols].copy()
    after = downholeData.shape[0]
    print(str(before-after)+' records removed without location information.')
    print(str(downholeData.shape[0])+' wells remain from '+str(downholeData['API_NUMBER'].unique().shape[0])+' located wells in study area.')
    return downholeData

#Function to remove headerData without surface topography information
##THIS ASSUMES AND SHOULD ONLY BE RUN AFTER ALL DESIRED SURFACE TOPO DATASETS HAVE BEEN MERGED/ADDED
def removenotopo(df, printouts=False):
    before = df.shape[0]
    
    df['ELEV_FT'].replace('', np.nan,inplace=True)
    df.dropna(subset=['ELEV_FT'], inplace=True)
    
    after = df.shape[0]

    if printouts:
        print("Number of rows before dropping those without surface elevation information: "+str(before))
        print("Number of rows after dropping those without surface elevation information: "+str(after))
        print('Well records deleted: '+str(before-after))
    return df

#This function drops all records in the downholedata with no depth information (either top or bottom depth of well interval)
def dropnodepth(df, printouts=False):
    #Replace empty cells in top and bottom columns with nan
    df['TOP'] = df['TOP'].replace('', np.nan)
    df['BOTTOM'] = df['BOTTOM'].replace('', np.nan)
    
    #Calculate number of rows before dropping
    before = df.shape[0]

    #Drop records without depth information
    df = df.dropna(subset=['TOP'])
    df = df.dropna(subset=['BOTTOM'])
    df.reset_index(inplace=True, drop=True) #Reset index
  
    if printouts:
        print("Number of rows before dropping those without record depth information: "+str(before))
        print("Number of rows after dropping those without record depth information: "+str(df.shape[0]))
        print('Number of well records without formation information deleted: '+str(before-df.shape[0])) 
    return df

#This function drops all records in downholeData with bad depth information (where the bottom of a record is nearer to the surface than the top)
def dropbaddepth(df, printouts=False):
    df['THICKNESS'] = df['BOTTOM'] - df['TOP'] #Calculate thickness
    before = df.shape[0] #Calculate number of rows before dropping
    df = df[(df['THICKNESS'] >= 0)] #Only include rows where thickness is positive (bottom is deeper than top)
    df.reset_index(inplace=True, drop=True) #Reset index

    if printouts:
        print("Number of rows before dropping those with obviously bad depth information: "+str(before))
        print("Number of rows after dropping those with obviously bad depth information: "+str(df.shape[0]))
        print('Well records deleted: '+str(before-df.shape[0]))
    return df

#This function drops all records in downholeData with no formation in formation in the description field
def dropnoformation(df, printouts=False):
    #Replace empty cells in formation column with nans
    df['FORMATION'] = df['FORMATION'].replace('', np.nan) 
    before = df.shape[0] #Calculate number of rows before dropping

    #Drop records without FORMATION information
    df = df.dropna(subset=['FORMATION'])
    df.reset_index(inplace=True, drop=True) #Reset index

    if printouts:
        print("Number of rows before dropping those without FORMATION information: "+str(before))
        print("Number of rows after dropping those without FORMATION information: "+str(df.shape[0]))
        print('Well records deleted: '+str(before-df.shape[0]))
        
    return df
