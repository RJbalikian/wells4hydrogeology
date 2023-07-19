# Wells 4 Hydrogeology (w4h) 

The w4h Python package is a package designed jointly by the Illinois State Geological Survey and Illinois State Water Survey.

It is designed to read in geology data from wells and create a layered, gridded hydrogeologic model of a study region, all within a python environment, automating and performing tasks often carried out in a dedicated GIS software.

The w4h package contains all the functions needed for getting N layers of a hydrogeology grid. 
Though the source code is split into separate modules, all functions are designed to be accessed directly from the w4h module (example: w4h.read_study_area())

The w4h module is designed to be flexible and customizable, allowing various kinds of data to be read in, with many different kinds of initial structures.

Using well descriptions from these database tables, the code contained here extracts, manipulates, and organizes the data to be used for hydrogeologic modeling. The scripts here can be used for specific regions of interests/study areas within the state, or for the state as a whole.

# API Documentation
API Documentation <a href="https://rjbalikian.github.io/wells4hydrogeology/main.html">here</a>

# Dependencies
## The w4h module has the following dependencies:
- numpy
- geopandas (and therefore pandas)
- rioxarray (and therefore xarray)
- matplotlib
- scipy
- owslib

# Inputs
## Required Inputs
Required inputs include are shown <a href="https://github.com/RJbalikian/wells4hydrogeology/">here</a>wiki

# Organization
## Modules
The package is organized by module, but all functions can be accessed directly using w4h.function_name() as well.
- core: general utility functions used throughout
- classify: functions for classifying the data
- clean: functions for cleaning the data
- export: functions for exporting the data, both as tables and rasters
- layers: functinos for generating layer(ed) models
- mapping: functions for mapping or performing geospatial analysis
- read: functions for reading in various files

# Included resources
The w4h package "ships" with some definition files that are read in/used by the scripts
- ISGS_HEADER_yyyy-mm-dd.TXT: tabular data exported from ISGS oracle database containing "header" information (i.e., metadata) about all the wells
- ISGS_DOWNHOLE_DATA_yyyy-mm-dd.TXT: tabular data exported from ISGS oracle database containing geologic information about wells in the state
- xyzData_yyyy-mm-dd.csv: most recent update of statewide wells with API, Latitude, Longitude, and surface elevation extracted from statewide lidar topography

# Intended workflow
Diagram for workflow available <a href="https://github.com/RJbalikian/wells4hydrogeology/wiki/Intended-Workflow-%5BDRAFT%5D">here</a>.