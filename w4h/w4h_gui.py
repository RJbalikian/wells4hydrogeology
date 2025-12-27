import copy
import inspect
import io
import pathlib
import sys
import tempfile

import geopandas as gpd
import matplotlib.pyplot as plt
from owslib.wms import WebMapService
import pandas as pd
import pyproj
import rioxarray as rxr
import streamlit as st
import xarray as xr

try:
    import w4h
    #from w4h import sprit_hvsr
    #from w4h import sprit_plot
except Exception:
    parent_dir = pathlib.Path(__file__).absolute().parent.parent.as_posix()
    sys.path.insert(0, parent_dir)
    import w4h


CRS_LIST = pyproj.database.query_crs_info()
CRS_STR_LIST = [f"{crs.name} ({crs.auth_name}:{crs.code})" for crs in CRS_LIST]
CRS_DICT = {f"{crs.name} ({crs.auth_name}:{crs.code})": crs for crs in CRS_LIST}
IL_LIDAR_URL=r"https://data.isgs.illinois.edu/arcgis/services/Elevation/IL_Statewide_Lidar_DEM_WGS/ImageServer/WMSServer?request=GetCapabilities&service=WMS"
GMRT_BASE_URL = r"https://www.gmrt.org:443/services/GridServer?minlongitude&maxlongitude%2C%20&minlatitude&maxlatitude&format=geotiff&resolution=default&layer=topo"

DEFAULT_POINTS_CRS = "WGS 84 (EPSG:4326)"
DEFAULT_POINTS_CRS_INDEX = CRS_STR_LIST.index(DEFAULT_POINTS_CRS)

RESOURCES = w4h.get_resources(scope='local')

VERBOSE = False

def get_defaults():
    w4h_funs = [w4h.file_setup, w4h.read_raw_csv, w4h.define_dtypes, w4h.merge_metadata, w4h.coords2geometry,
                w4h.read_study_area, w4h.clip_gdf2study_area, w4h.read_grid, w4h.add_control_points,
                w4h.remove_nonlocated, w4h.remove_no_topo, w4h.remove_no_depth, w4h.remove_bad_depth, w4h.remove_no_description,
                w4h.get_search_terms, w4h.read_dictionary_terms, w4h.specific_define, 
                w4h.split_defined, w4h.start_define, w4h.wildcard_define, w4h.remerge_data, w4h.fill_unclassified,
                w4h.read_lithologies, w4h.merge_lithologies, 
                w4h.align_rasters, w4h.get_drift_thick, w4h.sample_raster_points, w4h.get_layer_depths, w4h.layer_target_thick,
                w4h.layer_interp, w4h.export_grids]
    
    fullParamList = []
    fullDefaultList = []
    paramDefaultDict = {}
    for func in w4h_funs:
        parameters = inspect.signature(func).parameters
        defaults = [param.default for param in list(zip(*parameters.items()))[1]]
        parameters = list(zip(*parameters.items()))[0]
        fullParamList.extend(parameters)
        fullDefaultList.extend(defaults)

        for i, p in enumerate(parameters):
            paramDefaultDict[p] = defaults[i]
    return paramDefaultDict

st.session_state.param_defaults = param_defs = get_defaults()


def on_point_file_upload():
    file = st.session_state.point_file_ul
    temp_dir = tempfile.mkdtemp()
    path = pathlib.Path(temp_dir).joinpath(file.name)
    with open(path, "wb") as f:
            f.write(file.getvalue())
    if VERBOSE:
        print(file.name)
    st.session_state.well_data = path.as_posix()


def demo_checked():
    doDemo = st.session_state.demo_check


    well_data_file = RESOURCES['well_data']
    if doDemo:
        # Well data
        st.session_state.well_data = pathlib.Path(well_data_file).as_posix()

        # Surface raster
        st.session_state.surf_raster_type = 'Web Service'
        st.session_state.pre_service = 'ISGS Lidar'
        st.session_state.surf_raster_source = IL_LIDAR_URL

        # Bedrock raster
        st.session_state.lower_rast_TEXT = RESOURCES['bedrock_elev'].as_posix()

        # Model grid
        st.session_state.model_type = 'Raster Upload'
        st.session_state.model_grid_TEXT = RESOURCES['model_grid'].as_posix()

        # Study area
        st.session_state.study_area_TEXT = RESOURCES['study_area'].as_posix()


def surf_raster_type_changed():
    if st.session_state.surf_raster_type != "File":
        surf_raster_source_uploaded()
        service_changed()
    
    else: # st.session_state.surf_raster_type == "File":
        surf_raster_source_uploaded()
        service_changed()

        # Save service state
        if hasattr(st.session_state, 'pre_service'):
            st.session_state.preCServ = st.session_state.pre_service

        if hasattr(st.session_state, 'surf_rast_ul_obj_name'):
            st.session_state.surf_raster_source = st.session_state.surf_rast_ul_obj_name


def service_changed():
    if hasattr(st.session_state, 'preCServ') and not hasattr(st.session_state, 'pre_service'):
        srserve = st.session_state.preCServ
        if 'ISGS' in srserve:
            st.session_state.surf_raster_source = IL_LIDAR_URL
            st.session_state.surf_raster_CRS = "WGS 84 / Pseudo-Mercator (EPSG:3857)"
        elif "GMRT" in srserve:
            st.session_state.surf_raster_source = GMRT_BASE_URL
            st.session_state.surf_raster_CRS = "WGS 84 (EPSG:4326)"
    elif hasattr(st.session_state, 'pre_service'):
        st.session_state.preCServ = st.session_state.pre_service
        srserve = st.session_state.preCServ
        if 'ISGS' in srserve:
            st.session_state.surf_raster_source = IL_LIDAR_URL
            st.session_state.surf_raster_CRS = "WGS 84 / Pseudo-Mercator (EPSG:3857)"
        elif "GMRT" in srserve:
            st.session_state.surf_raster_source = GMRT_BASE_URL
            st.session_state.surf_raster_CRS = "WGS 84 (EPSG:4326)"
    
    on_surf_raster_source_change()


def surf_raster_source_uploaded():
    if hasattr(st.session_state, 'surf_rast_ul'):
        st.session_state.surf_rast_ul_obj = copy.deepcopy(st.session_state.surf_rast_ul)
    on_surf_raster_source_change()
 
    
def on_surf_raster_source_change():
    rType = st.session_state.surf_raster_type
    rsource = st.session_state.surf_raster_source

    #st.write("ON SURF RAST SOURCE CHAGNE")
    #st.write(st.session_state.surf_raster_type)

    if rType == "File":
        print("FILE")
        if not hasattr(st.session_state, 'surf_rast_ul') and  hasattr(st.session_state, 'surf_rast_ul_name'):
            st.session_state.surf_raster_source = pathlib.Path(st.session_state.surf_rast_ul_name).as_posix()

        if hasattr(st.session_state, 'surf_rast_ul') and st.session_state.surf_rast_ul is not None:
            st.write(pathlib.Path(st.session_state.surf_rast_ul.name).as_posix())
            st.session_state.surf_raster_source = pathlib.Path(st.session_state.surf_rast_ul.name).as_posix()
        elif hasattr(st.session_state, 'surf_rast_ul_name'):
            st.session_state.surf_raster_source = pathlib.Path(st.session_state.surf_rast_ul_name).as_posix()
        else:
            st.session_state.surf_raster_source = 'None'
    else:
        if hasattr(st.session_state, 'surf_rast_ul') and st.session_state.surf_rast_ul is not None:
            st.session_state.surf_rast_ul_name = st.session_state.surf_rast_ul.name

        if hasattr(st.session_state, 'pre_service'):
            srserve = st.session_state.pre_service
            if 'ISGS' in srserve:
                st.session_state.pre_service = 'ISGS Lidar'
                st.session_state.surf_raster_CRS = "WGS 84 / Pseudo-Mercator (EPSG:3857)"
                st.session_state.surf_raster_source = IL_LIDAR_URL
            elif "GMRT" in srserve:
                st.session_state.pre_service = 'GMRT'
                st.session_state.surf_raster_source = GMRT_BASE_URL
                st.session_state.surf_raster_CRS = "WGS 84 (EPSG:4326)"


def w4hrun():
    stss = st.session_state
    st.toast('Processing')

    specified_params_dict = {'well_data':'well_data',
                        'surf_elev_grid':'surf_raster_source', 
                        'bedrock_elev_grid':'lower_rast_TEXT', 
                        'model_grid':'model_grid_TEXT',
                        'study_area':'study_area_TEXT'}
    # Code to read in well_data file
    # FIX THIS

    #if st.session_state.demo_check:
    #    RESOURCES = w4h.get_resources(scope='local')
    #    st.session_state.well_data = RESOURCES['well_data'].as_posix()


    surf_elev = None
    if stss.surf_raster_type == 'File':
        pass
        # Code to read in file
        stss.surf_elev_grid = None
    else:
        # Code to read in service
        if 'ISGS' in st.session_state.pre_service:
            elevation_source = st.session_state.surf_raster_source
            wms = WebMapService(elevation_source)

            bboxFile = gpd.read_file(RESOURCES['study_area'])
            bboxFile = bboxFile.to_crs('EPSG:4326')
            bboxFile = bboxFile.bounds
            bbox = bboxFile.values.tolist()[0]

            img = wms.getmap(
                layers=['IL_Statewide_Lidar_DEM_WGS:None'],
                srs='EPSG:4326',
                bbox=bbox,
                size=(256, 256),
                format='image/tiff',
                transparent=True
                )

            bio = io.BytesIO(img.read())
            elevData_rxr = rxr.open_rasterio(bio)
            elevData_rxr = elevData_rxr.rio.write_crs('EPSG:4326')
            #elevData_rxr = elevData_rxr.isel(band=0)

            output_crs = "EPSG:4326"
            surf_elev = elevData_rxr.rio.reproject(output_crs)
        else:
            surf_elev = None

    # Code to read in lower surface elevation
    # Fix this to get it working with actual file, not just text
    br_elev = None
    if '.tif' in str(stss.lower_rast_TEXT):
        stss.lower_rast_TEXT = pathlib.Path(stss.lower_rast_TEXT).with_suffix('.TIF').as_posix()
        br_elev = rxr.open_rasterio(stss.lower_rast_TEXT).rio.reproject('EPSG:4326')

    if '.tif' in str(stss.model_grid_TEXT):
        stss.model_grid_TEXT = pathlib.Path(stss.model_grid_TEXT).with_suffix('.TIF').as_posix()


    # Get model grid
    if 'node' in str(stss.model_type).lower():
        # FIX code to generate model grid from node spacing or number
        pass
    elif stss.model_type == 'Surface Elevation':
        # CODE TO copy SURFACE ELEVATION grid
        pass
    elif stss.model_type == 'Lower Surface':
        # CODE TO copy SURFACE ELEVATION grid
        pass
    else:
        # Code to read raster
        # FIX THIS
        stss.model_grid = rxr.open_rasterio(stss.model_grid_TEXT)
    
    # Add non-default params to kwargs
    w4hrun_kwargs = {}
    for paramName, defaultVal in stss.param_defaults.items():
        if paramName not in specified_params_dict:
            if hasattr(st.session_state, paramName):
                if stss[paramName] != defaultVal:
                    w4hrun_kwargs[paramName] = stss[paramName]
    
    # (Eventually) Code to print out what is not default value

    # Specify unmodifiable kwargs
    if surf_elev is not None:
        del specified_params_dict['surf_elev_grid']
        w4hrun_kwargs['surf_elev_grid'] = surf_elev
        #st.write(surf_elev['spatial_ref'])
    if br_elev is not None:
        del specified_params_dict['bedrock_elev_grid']
        w4hrun_kwargs['bedrock_elev_grid'] = br_elev

    for kw, stsskw in specified_params_dict.items():
        w4hrun_kwargs[kw] = stss[stsskw]
    w4hrun_kwargs['verbose'] = True

    if 'output_crs' in w4hrun_kwargs.keys():
        w4hrun_kwargs['output_crs'] = 'EPSG'+w4hrun_kwargs['output_crs'].split('EPSG')[1][:-1]

    # Run the processing
    with st.spinner('Processing data', show_time=True):
        #st.write(w4hrun_kwargs)
        with st.status('Running'):
            try:
                st.session_state.results = w4h.run(**w4hrun_kwargs)
            except Exception as e:
                st.error(e)
                print(e)
                st.session_state.results = None
    st.success('Processing complete')

    resDS = st.session_state.results[1]
    if resDS.coords['x'][0] > resDS.coords['x'][-1]:
        resDS = resDS.sel(x=resDS.x[::-1])

    if resDS.coords['y'][0] > resDS.coords['y'][-1]:
        resDS = resDS.sel(y=resDS.y[::-1])
    
    #st.session_state.results[1] = resDS
    #st.write(st.session_state.results)
    fig, ax = plt.subplots(figsize=(10, 10))
    #st.session_state.results[1]['Model_Layers'][0].plot(ax=ax)

    if stss['study_area_TEXT'] is not None:
        studyArea = gpd.read_file(stss['study_area_TEXT'])
        studyArea = studyArea.to_crs(w4hrun_kwargs['output_crs'])
        plotDA = resDS['Model_Layers'].sum(axis=0)

        plotDARP = plotDA.rio.reproject(w4hrun_kwargs['output_crs'])
        plotDARP = plotDARP.rio.clip(studyArea.geometry, studyArea.crs)

        plotDARP.plot(ax=ax)
        studyArea.plot(ax=ax, facecolor="#00000000", edgecolor='k')

    st.pyplot(fig)


def main():

    st.set_page_config(page_title='W4H WebApp',
                       page_icon=":material/globe_book:",
                       layout='wide',
                       initial_sidebar_state='expanded',
                       menu_items={"Get Help":'https://github.com/RJbalikian/wells4hydrogeology/wiki',
                                    "Report a bug":"https://github.com/RJbalikian/wells4hydrogeology/issues",
                                    'About':"W4H is a collaboration between the Illinois State Water Survey and the Illinois State Geological Survey"})
    with st.sidebar:
        sampleCol, headerCol  = st.columns([0.7, 0.3], vertical_alignment='top')
        with headerCol.container(horizontal_alignment='right'):
            st.button('Run Analysis', type='primary', on_click=w4hrun, key='run_button')
        with sampleCol.container(horizontal_alignment='right'):
            st.checkbox('Demo run', disabled=False, value=False, key='demo_check', 
                        on_change=demo_checked)
        st.header("Specify Input Data", divider='rainbow')
        

        with st.expander("Well Data", expanded=True):
            wdval = None
            if hasattr(st.session_state, 'point_file_ul') and st.session_state.point_file_ul is not None:
                wdval = st.session_state.point_file_ul.name
            st.text_input(label="Well data file", value=wdval, key='well_data')
            st.file_uploader(label='Upload Point File', key='point_file_ul',
                            on_change=on_point_file_upload)

        with st.expander("Raster Data", expanded=True):
            st.header('Raster Data')
            surfTab, brtab = st.tabs(["Surface Raster", "Lower Raster"])
            surfTab.segmented_control(label='Select Raster Type', options=['File', 'Web Service'], 
                                      key='surf_raster_type', default='Web Service',
                                      on_change=surf_raster_type_changed)

            if st.session_state.surf_raster_type == 'File':
                srval = None
                if hasattr(st.session_state, 'surf_rast_ul') and st.session_state.surf_rast_ul is not None:
                    srval = st.session_state.surf_rast_ul.name
                surfTab.text_input(label="Surface Raster file", value=srval, key='surf_raster_source',
                                    on_change=on_surf_raster_source_change)
                surfTab.file_uploader(label='Upload Surface Elevation Raster', key='surf_rast_ul',
                                        on_change=surf_raster_source_uploaded)

            else:
                srsInd = 0
                if hasattr(st.session_state, 'preCServ') and not hasattr(st.session_state, 'pre_service'):
                    srserve = st.session_state.preCServ
                    if 'ISGS' in srserve:
                        srsInd = 1
                    elif "Custom" in srserve:
                        srsInd = 2

                surfTab.radio(label='Surface Raster Services', 
                            options=['GMRT', 'ISGS Lidar', "Custom"], index=srsInd,
                            horizontal=True, 
                            key='pre_service', on_change=service_changed)
                servTextDisabled = False
                if st.session_state.pre_service != 'Custom':
                    servTextDisabled = True
                
                srval = None
                if "ISGS" in st.session_state.pre_service:
                    srval = IL_LIDAR_URL
                elif "GMRT" in st.session_state.pre_service:
                    srval = GMRT_BASE_URL

                surfTab.text_input(label='Specify Surface Raster Service URL (Currently only WMS supported)',
                                   value=srval,
                                   key='surf_raster_source')

                if hasattr(st.session_state, "surf_service_url"):
                    if st.session_state.pre_service == "GMRT":
                        st.session_state.surf_service_url = GMRT_BASE_URL
                    elif st.session_state.pre_service == "ISGS Lidar":
                        st.session_state.surf_service_url = IL_LIDAR_URL



            specSurfRastCol, surfRastCRSCol = surfTab.columns([0.3, 0.7])
            specSurfRastCol.toggle('Specify Surface Raster CRS', key='specify_surfrast_crs')

            surfRastCRSDisabled =  not st.session_state.specify_surfrast_crs
            surfRastCRSCol.selectbox("Surface Raster CRS", disabled = surfRastCRSDisabled,
                              options=CRS_STR_LIST, index=DEFAULT_POINTS_CRS_INDEX, 
                              key='surf_raster_CRS')

            # Bedrock raster
            brval = None
            if hasattr(st.session_state, 'lower_rast_UL') and st.session_state.lower_rast_UL is not None:
                brval = st.session_state.lower_rast_UL.name
            if hasattr(st.session_state, 'lower_rast_TEXT'):
                st.session_state.lower_rast_TEXT = pathlib.Path(st.session_state.lower_rast_TEXT).as_posix()
                brval = st.session_state.lower_rast_TEXT
            brtab.text_input(label="Lower Raster File", value=brval, key='lower_rast_TEXT')
            brtab.file_uploader(label='Upload Lower Elevation Raster', key='lower_rast_UL')

        with st.expander("Extent and Resolution"):
            st.selectbox('Output CRS', options=CRS_STR_LIST, 
                         index=DEFAULT_POINTS_CRS_INDEX,
                         key='output_crs')

            saval = None
            if hasattr(st.session_state, 'study_area_ul') and st.session_state.lower_rast_UL is not None:
                saval = st.session_state.study_area_ul.name
            st.text_input(label="Study Area File", value=saval, key='study_area_TEXT')            
            st.file_uploader('Upload Study Area File', key='study_area_ul')
            st.header('Model Grid')
            st.selectbox("Model Grid Source", 
                         options=["Lower Surface", "Surface Elevation",
                                  "Raster Upload", "# Nodes", "Node Size"],
                         index=0, key='model_type')
            if 'node' not in st.session_state.model_type.lower():
                if st.session_state.model_type == 'Raster Upload':
                    mgval = None
                    if hasattr(st.session_state, 'model_grid_UL') and st.session_state.model_grid_UL is not None:
                        mgval = st.session_state.model_grid_UL.name
                    st.text_input(label='Model Grid File', value=mgval, key='model_grid_TEXT')
                    st.file_uploader(label='ModelGrid', key='model_grid_UL')
            else:
                xnodecol, ynodecol = st.columns([0.5, 0.5])
                xnodeLabel = 'No. X Nodes'
                ynodeLabel = 'No. Y Nodes'
                if st.session_state.model_type == 'Node Size':
                    xnodeLabel = 'X Node Size'
                    ynodeLabel = 'Y Node Size'

                xnodecol.number_input(xnodeLabel, min_value=5, max_value=5000,
                                      step=1, value=100, key='x_node_step')
                ynodecol.number_input(ynodeLabel, min_value=5, max_value=5000,
                                      step=1, value=100, key='y_node_step')
        
        with st.expander("Data Mapping", expanded=True):
            widcol, descol = st.columns([0.5, 0.5])
            widcol.selectbox(label="Well ID Column", options=['API_NUMBER'], 
                            help="Select name of column from point file containing unique well ID's", key='well_id_col')
            descol.selectbox(label="Description Column", options=['FORMATION'], 
                            help="Select name of column from point file containing lithologic descriptions", key='description_col')

            xcol, ycol, zcol = st.columns([0.3, 0.3, 0.3])
            xcol.selectbox(label="X Coord Column", options=['LONGITUDE'])
            ycol.selectbox(label="Y Coord Column", options=['LATITUDE'])
            zcol.selectbox(label="Elevation Column", options=['SURFACE_ELEV'])

            tcol, bcol = st.columns([0.5, 0.5])
            tcol.selectbox(label="Top Column", options=['TOP', '2'])
            bcol.selectbox(label="Bottom Column", options=['BOTTOM', '2'])

        with st.expander("Additional Settings"):
            wellSetTab, surfElevSetTab, botSetTab = st.tabs(["WellData", "Surface Elevation", "Bottom Layer"])

if __name__ == "__main__":
    main()
