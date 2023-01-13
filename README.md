# Wells4hydrogeology
This repository uses data from multiple tables in the statewide well database of the Illinois State Geological Survey (ISGS). 

Using well descriptions from these database tables, the code contained here extracts, manipulates, and organizes the data to be used for hydrogeologic modeling. The scripts here can be used for specific regions of interests/study areas within the state, or for the state as a whole.

Required inputs include:
- ISGS_DOWNHOLE_DATA: A table in the ISGS database containing descriptions of well intervals for wells throughout the state.
- ISGS_HEADER: A table in the ISGS database containing the metadata for all the wells, including API number, well location, and in some cases elevation.
- XYZData (optional): A separate table containing updated location information for each well, particularly updated elevation data from Lidar products
- Surface elevation: raster data containing the surface elevation of the study area or state.
- Bedrock elevation: raster data containing the bedrock elevation of the study area or state.
- Model grid: raster data whose resolution and cell locations align with that of the hydrogeologic model (i.e., in MODFLOW)
- Lithology data: "definitions" to convert raw, manual well descriptions to broad lithological categories and then to a target lithology (e.g., coarse sediment).

## master_notebook contains an interactive jupyter notebook with all the steps for running the main body of the script

## lib folder contains all scripts with functions used
- readData: functions for reading in various files
- mapping: functions for mapping or performing geospatial analysis
- cleanData: functions for cleaning the data
- classify: functions for classifying the data

## res folder contains all the files that are read in/used by the scripts
- ISGS_HEADER_yyyy-mm-dd.TXT: tabular data exported from ISGS oracle database containing "header" information (i.e., metadata) about all the wells
- ISGS_DOWNHOLE_DATA_yyyy-mm-dd.TXT: tabular data exported from ISGS oracle database containing geologic information about wells in the state
- xyzData_yyyy-mm-dd.csv: most recent update of statewide wells with API, Latitude, Longitude, and surface elevation extracted from statewide lidar topography

## doc folder contains documentation 

## out folder contains outputs from script runs
