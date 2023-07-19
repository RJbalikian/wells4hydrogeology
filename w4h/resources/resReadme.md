# Resources in this folder

## Datatype dictionaries
**Files with dictionaries for assigning proper numpy datatypes to dataframes from data that is read in**
- *donwholeDataTypes.txt*: contains a dictionary of datatypes for the ISGS_donwholeData Oracle table when it is read in as a pandas dataframe
- *headerDataTypes.txt*: contains a dictionary of datatypes for the ISGS_header Oracle table when it is read in as a pandas dataframe
- *xyzDataTypes.txt*: contains a dictionary of datatypes for the xyzData csv file when it is read in as a pandas dataframe (this is currently generated outside this code)
- *isws_crs.txt*: contains a dictionary with all the relevant information for use with pyproj to create a well-known text (WKT) string of the CRS

## Lithology Interpretations
**Files to convert the automated geologic classifications into lithology for the various kinds of interpretation used in the output for hydrogeologic parameterization**
*Lithology_Interp_...*:
- *...Clay.csv*: Contains "dummy variables" for the interpretations that should be classified as clay (or not)
- *...FineCoarse.csv*: Contains "dummy variables" for the interpretations that should be classified as coarse (or not) (this is for broad classification)
- *...Gravel.csv*: Contains "dummy variables" for the interpretations that should be classified as gravel (or not)
- *...Sand.csv*: Contains "dummy variables" for the interpretations that should be classified as sand (or not)
- *...Silt.csv*: Contains "dummy variables" for the interpretations that should be classified as silt (or not)

## Search terms
**Files containing search terms used to classify descriptions made by geologists/in the field to a subset of geologic materials**
*SearchTerms-...*:
- *...Specific_...csv*: Contains all of the geologic descriptions from the downholeData table that have been specificallly classified by a geologist. 
- *...Start.csv*: Contains geologic descriptions that are less common, but which whose lithologies can be interprted from the first substring of the description
    - (For example, if a description begins with "limestone, ..." it is highly likely that, no matter how the rest of the description modifies this, this interval can be classified as bedrock)
