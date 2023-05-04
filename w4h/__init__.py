#__init__.py
"""The Wells 4 Hydrogeology (w4h) Python package is a package designed jointly by the Illinois State Geological Survey and Illinois State Water Survey.

It is designed to read in geology data from wells and create a layered, gridded hydrogeologic model of a study region, all within a python environment, automating and performing tasks often carried out in a dedicated GIS software.

The w4h package contains all the functions needed for getting N layers of a hydrogeology grid. 
Though the source code is split into separate modules, all functions are designed to be accessed directly from the w4h module (example: w4h.read_study_area())

The w4h module is designed to be flexible and customizable, allowing various kinds of data to be read in, with many different kinds of initial structures.

The w4h module has the following dependencies:
- numpy
- pandas
- rioxarray (and therefore xarray)
- geopandas
- matplotlib
- scipy
- owslib

"""

#from w4h import classify, clean, export, layers, mapping, read

from w4h.utilities import(logger)

from w4h.classify import (specific_define, 
                          split_defined, 
                          start_define, 
                          remerge_data, 
                          depth_define, 
                          export_undefined, 
                          fill_unclassified, 
                          merge_lithologies, 
                          get_unique_wells,
                          sort_dataframe)

from w4h.clean import (remove_nonlocated, 
                       remove_no_topo, 
                       drop_no_depth, 
                       drop_bad_depth, 
                       drop_no_formation)

from w4h.export import (export_dataframe,
                        export_grids)

from w4h.layers import (get_layer_depths,
                        merge_tables, 
                        layer_target_thick, 
                        layer_interp,
                        combine_dataset)

from w4h.mapping import (read_study_area, 
                         coords2geometry, 
                         clip_gdf2study_area, 
                         sample_raster_points, 
                         xyz_metadata_merge, 
                         read_wms,
                         read_wcs, 
                         grid2study_area,
                         read_model_grid,
                         read_grid,
                         align_rasters,
                         get_drift_thick)

from w4h.read import (get_current_date,
                      get_most_recent,
                      file_setup,
                      read_raw_txt,
                      read_xyz,
                      read_dict,
                      define_dtypes,
                      get_search_terms,
                      read_dictionary_terms,
                      read_lithologies)


__all__=(
        'specific_define', 
        'split_defined', 
        'start_define', 
        'remerge_data', 
        'depth_define', 
        'export_undefined', 
        'fill_unclassified', 
        'merge_lithologies', 
        'get_unique_wells',
        'sort_dataframe',
         'remove_nonlocated', 
         'remove_no_topo', 
         'drop_no_depth', 
         'drop_bad_depth', 
         'drop_no_formation',
        'export_dataframe',
        'export_grids',
         'get_layer_depths',
         'merge_tables', 
         'layer_target_thick', 
         'layer_interp',
         'combine_dataset',
        'read_study_area', 
        'coords2geometry', 
        'clip_gdf2study_area', 
        'sample_raster_points', 
        'xyz_metadata_merge', 
        'read_wms',
        'read_wcs', 
        'grid2study_area',
        'read_model_grid',
        'read_grid',
        'align_rasters',
        'get_drift_thick',
         'get_current_date',
         'get_most_recent',
         'file_setup',
         'read_raw_txt',
         'read_xyz',
         'read_dict',
         'define_dtypes',
         'get_search_terms',
         'read_dictionary_terms',
         'read_lithologies',
        'logger')

__author__='Riley Balikian, Joe Franke, Allan Jones, Mike Krasowski'