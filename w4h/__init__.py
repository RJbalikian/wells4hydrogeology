#__init__.py
"""This module contains all the functions needed for getting 9 layers of geology"""

from w4h import classify, clean, export, layers, mapping, read

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
                         coords2Geometry, 
                         clipHeader2StudyArea, 
                         sample_raster_points, 
                         addElevtoHeader, 
                         readWMS,
                         readWCS, 
                         clipGrid2StudyArea,
                         read_model_grid,
                         read_grid,
                         alignRasters,
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
        'classify', 
         'clean', 
         'export', 
         'layers',
         'mapping',
         'read',
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
        'coords2Geometry', 
        'clipHeader2StudyArea', 
        'sample_raster_points', 
        'addElevtoHeader', 
        'readWMS',
        'readWCS', 
        'clipGrid2StudyArea',
        'read_model_grid',
        'read_grid',
        'alignRasters',
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
         'read_lithologies')

__author__='Riley Balikian, Joe Franke, Allan Jones, Mike Krasowski'