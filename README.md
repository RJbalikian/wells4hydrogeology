# Wells4hydrogeology
Using well descriptions from database, create data to be used for groundwater hydrogeological model

## master_notebook contains an interactive jupyter notebook with all the steps for running the main body of the script

## lib folder contains all scripts with functions used
- setupFiles: functions for setting up filepaths, variables, etc.
- readData: functions for reading in various files
- cleanData: functions for cleaning the data
- classify: functions for classifying the data
- mapping: functions for mapping or performing geospatial analysis

## res folder contains all the files that are read in/used by the scripts
- ISGS_HEADER_yyyy-mm-dd.TXT: tabular data exported from ISGS oracle database containing "header" information (i.e., metadata) about all the wells
- ISGS_DOWNHOLE_DATA_yyyy-mm-dd.TXT: tabular data exported from ISGS oracle database containing geologic information about wells in the state
- xyzData_yyyy-mm-dd.csv: most recent update of statewide wells with API, Latitude, Longitude, and surface elevation extracted from statewide lidar topography

## doc folder contains documentation 

##out folder contains outputs from script runs
