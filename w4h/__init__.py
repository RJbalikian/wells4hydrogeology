#__init__.py
"""This module contains all the functions needed for getting 9 layers of geology"""

from w4h import classify, clean, export, layers, mapping, read


from w4h.classify import (specificDefine, 
                          splitDefined, 
                          startDefine, 
                          remergeData, 
                          depthDefine, 
                          export_toBeDefined, 
                          fillUnclassified, 
                          mergeLithologies, 
                          getUniqueWells)


from w4h.clean import (removeNonlocatedData, 
                       removenotopo, 
                       dropnodepth, 
                       dropbaddepth, 
                       dropnoformation)

from w4h.export import (exportDataframe)

from w4h.layers import (get_layer_depths,
                        merge_tables, 
                        layer_target_thick, 
                        layer_interp)

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

from w4h.read import (getCurrentDate,
                      findMostRecentFiles,
                      filesSetup,
                      readRawTxtData,
                      readXYZData,
                      get_resource_path,
                      read_dict,
                      defineDataTypes,
                      searchTermFilePaths,
                      readSearchTerms,
                      readLithologies)

__all__=(
        'classify', 
         'clean', 
         'export', 
         'layers',
         'mapping',
         'read',
        'specificDefine', 
        'splitDefined', 
        'startDefine', 
        'remergeData', 
        'depthDefine', 
        'export_toBeDefined', 
        'fillUnclassified', 
        'mergeLithologies', 
        'getUniqueWells',
         'removeNonlocatedData', 
         'removenotopo', 
         'dropnodepth', 
         'dropbaddepth', 
         'dropnoformation',
        'exportDataframe',
         'get_layer_depths',
         'merge_tables', 
         'layer_target_thick', 
         'layer_interp',
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
         'getCurrentDate',
         'findMostRecentFiles',
         'filesSetup',
         'readRawTxtData',
         'readXYZData',
         'get_resource_path',
         'read_dict',
         'defineDataTypes',
         'searchTermFilePaths',
         'readSearchTerms',
         'readLithologies')

__author__='Riley Balikian, Joe Franke'