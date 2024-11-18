# __init__.py
"""
This is the wells4hydrogeology package.
  
It contains the functions needed to convert raw well descriptions into usable (hydro)geologic data.

"""
__version__ = "0.0.25"


from w4h.core import (logger_function,
                      verbose_print,
                      run,
                      get_resources)

from w4h.classify import (specific_define,
                          split_defined,
                          start_define,
                          wildcard_define,
                          remerge_data,
                          depth_define,
                          export_undefined,
                          fill_unclassified,
                          merge_lithologies,
                          get_unique_wells,
                          sort_dataframe)

from w4h.clean import (remove_nonlocated,
                       remove_no_topo,
                       remove_no_depth,
                       remove_bad_depth,
                       remove_no_description)

from w4h.export import (export_dataframe,
                        export_grids)

from w4h.layers import (get_layer_depths,
                        merge_metadata,
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
                      read_raw_csv,
                      read_xyz,
                      read_dict,
                      define_dtypes,
                      get_search_terms,
                      read_dictionary_terms,
                      read_lithologies,
                      add_control_points)

from w4h.visualization import(
                            plot_cross_section
                             )

__all__ = ('logger_function', 'verbose_print', 'run', 'get_resources',
           # Classify
           'specific_define',
           'split_defined',
           'start_define',
           'wildcard_define',
           'remerge_data',
           'depth_define',
           'export_undefined',
           'fill_unclassified',
           'merge_lithologies',
           'get_unique_wells',
           'sort_dataframe',
           # Clean
           'remove_nonlocated',
           'remove_no_topo',
           'remove_no_depth',
           'remove_bad_depth',
           'remove_no_description',
           # Export
           'export_dataframe',
           'export_grids',
           # Layers
           'get_layer_depths',
           'merge_metadata',
           'layer_target_thick',
           'layer_interp',
           'combine_dataset',
           # Mapping
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
           # Read
           'get_current_date',
           'get_most_recent',
           'file_setup',
           'read_raw_csv',
           'read_xyz',
           'read_dict',
           'define_dtypes',
           'get_search_terms',
           'read_dictionary_terms',
           'read_lithologies',
           'add_control_points',
           # Visualization
           'plot_cross_section'
           )

# Update the w4h.run() help() return to actually be helpful
run.__doc__ = core._run_docstring()
__author__ = 'Riley Balikian, \
              Joe Franke, Mike Krasowski, Allan Jones, Daniel Abrams'
