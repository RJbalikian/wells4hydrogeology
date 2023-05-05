# Wells 4 Hydrogeology (w4h) 

The w4h Python package is a package designed jointly by the Illinois State Geological Survey and Illinois State Water Survey.

It is designed to read in geology data from wells and create a layered, gridded hydrogeologic model of a study region, all within a python environment, automating and performing tasks often carried out in a dedicated GIS software.

The w4h package contains all the functions needed for getting N layers of a hydrogeology grid. 
Though the source code is split into separate modules, all functions are designed to be accessed directly from the w4h module (example: w4h.read_study_area())

The w4h module is designed to be flexible and customizable, allowing various kinds of data to be read in, with many different kinds of initial structures.

Using well descriptions from these database tables, the code contained here extracts, manipulates, and organizes the data to be used for hydrogeologic modeling. The scripts here can be used for specific regions of interests/study areas within the state, or for the state as a whole.

# API Documentation
<a href="main.html">API Documentation</a>

# Dependencies
 The w4h module has the following dependencies:
- numpy
- pandas
- rioxarray (and therefore xarray)
- geopandas
- matplotlib
- scipy
- owslib

# Inputs
Required inputs include:
- ISGS_DOWNHOLE_DATA: A table in the ISGS database containing descriptions of well intervals for wells throughout the state.
- ISGS_HEADER: A table in the ISGS database containing the metadata for all the wells, including API number, well location, and in some cases elevation.
- XYZData (optional): A separate table containing updated location information for each well, particularly updated elevation data from Lidar products
- Surface elevation: raster data containing the surface elevation of the study area or state.
- Bedrock elevation: raster data containing the bedrock elevation of the study area or state.
- Model grid: raster data whose resolution and cell locations align with that of the hydrogeologic model (i.e., in MODFLOW)
- Lithology data: "definitions" to convert raw, manual well descriptions to broad lithological categories and then to a target lithology (e.g., coarse sediment).

## master_notebook contains an interactive jupyter notebook with all the steps for running the main body of the script

## w4h folder contains all scripts with functions used
- readData: functions for reading in various files
- mapping: functions for mapping or performing geospatial analysis
- cleanData: functions for cleaning the data
- classify: functions for classifying the data

## resources folder contains all the files that are read in/used by the scripts
- ISGS_HEADER_yyyy-mm-dd.TXT: tabular data exported from ISGS oracle database containing "header" information (i.e., metadata) about all the wells
- ISGS_DOWNHOLE_DATA_yyyy-mm-dd.TXT: tabular data exported from ISGS oracle database containing geologic information about wells in the state
- xyzData_yyyy-mm-dd.csv: most recent update of statewide wells with API, Latitude, Longitude, and surface elevation extracted from statewide lidar topography

## Intended workflow 

```mermaid
flowchart TD
    subgraph setup["Set up Files"]
        direction RL

        A0((file_setup))

        A1[(db_dir)]
        A2{{log_dir}}
        A3{{verbose}}
        A4{{log}}

        A1-->A0
        A2-.->A0
        A3-.->A0
        A4-.->A0
    end

    subgraph read["Read Data"]
        direction TB
        subgraph readraw["read_raw_txt()"]
            direction RL
            B0(("read_raw_txt()"))
            B1[(filepath)]
            B2{{"pd.read_csv(**kwargs)"}}

            B1-->B0
            B2-.->B0
        end

        subgraph defineDtypes["define_dtypes()"]
            direction RL
            C0(("define_datatypes()"))
            C1[(df)]
            C2{{"dtypes"}}

            C1-->C0
            C2-->C0
        end
        readraw-->defineDtypes
    end

    subgraph readSA["Read Study Area"]
        direction RL
        D0(("read_study_area()"))
        D1[(studyareapath)]
        
        D1---D0

    end

    subgraph Read_Grids["Read Grids"]
        direction RL
        E0(("read_grids()"))
        E1[(datapath)]
        E2{{grid_type}}
        E3[(study_area)]
        E4{{clip_to_studyarea}}
        E5{{use_service}}

        E1-->E0
        E2-.->E0
        E3-.->E0
        E4-.->E0
        E5-.->E0

    end

    subgraph clipdatatoSA["Clip Data to Study Area"]
        direction TB
        subgraph xyzmetaMerge
            direction RL
            F0(("xyz_metadata_merge()"))
            F1((xyz))
            F2((metadata))

            F1-->F0
            F2-.->F0
        end

        subgraph coords2geometry
            direction RL
            G0(("coords2geometry()"))
            G1[(df)]
            G2{{xcol}}
            G3{{ycol}}
            G4{{zcol}}
            G5{{crs}}

            G1-->G0
            G2-.->G0
            G3-.->G0
            G4-.->G0
            G5-.->G0

        end

        subgraph clipdata[Clip Data]
            direction RL
            H0(("clip_gdf2study_area()"))  
            H1[("study_area")]
            H2[("gdf")]
            H3{{"gdf_crs"}}

            H1-->H0
            H2-.->H0
            H3-.->H0
        end

        xyzmetaMerge --> coords2geometry  
        coords2geometry --> clipdata  
    end

    subgraph clean["Clean Data"]
        direction TB

        subgraph nonlocated
            direction RL
            I0(("remove_nonlocated()"))
            I1[(df)]
            I2[(metadata_df)]

            I1-->I0
            I2-->I0
        end

        subgraph notopo
            direction RL
            J0(("remove_no_topo()"))
            J1[(df)]

            J1-->J0
        end

        subgraph nodepth
            direction RL
            K0(("drop_no_depth()"))
            K1[(df)]

            K1-->K0
        end

        subgraph baddepth
            direction RL
            L0(("drop_bad_depth()"))
            L1[(df)]

            L1-->L0
        end

        subgraph noformation
            direction RL
            M0(("drop_no_formation()"))
            M1[(df)]

            M1-->M0
        end

        subgraph merge["MERGE:NEEDED?"]
            direction RL
            N0((merge))
        end

        nonlocated--> notopo
        notopo-->nodepth
        nodepth-->baddepth
        baddepth-->noformation
        noformation-->merge
    end

    subgraph classify["Classify Well Descriptions"]
        direction TB

        subgraph getterms
            direction RL
            O0(("get_search_terms()"))
            O1[("spec_dir()")]
            O2{{"spec_glob_pattern()"}}
            O3[("start_dir()")]
            O4{{"start_glob_pattern()"}}

            O1-->O0
            O2-->O0
            O3-.->O0
            O4-->O0
        end

        subgraph readterms
            direction RL
            P0(("read_dictionary_terms()"))
            P1{{dict_file}}

            P1-->P0
        end

        readtermsIter[["Iterate read_terms() as needed"]]

        subgraph define[Define]
            direction LR

            subgraph specdefine['Define Exact Matches]
                direction RL
                Q0(("specific_define()"))
                Q1[(df)]
                Q2[(terms_df)]

                Q1-->Q0
                Q2-->Q0
            end

            subgraph startdefine["Define Substring Matches"]
                direction RL
                R0(("start_define()"))
                R1[(df)]
                R2[(terms_df)]

                R1-->R0
                R2-->R0                
            end

            subgraph depthdefine["Define by Depth"]
                direction RL
                S0(("depth_define()"))
                S1[(df)]
                S2[(thresh)]

                S1-->S0
                S2-->S0                
            end

            specdefine~~~startdefine
            startdefine~~~depthdefine
        end

        subgraph fillunclass
            direction RL
            T0(("fill_unclassified()"))
            T1[(df)]

            T1-->T0
        end
        
        subgraph readlith
            direction RL
            U0(("read_lithologies()"))
        end

        subgraph mergelith
            direction RL
            V0(("merge_lithologies()"))
            V1[(df)]
            V2[(targinterps_df)]

            V1-->V0
            V2-->V0
        end

        getterms-.->readterms
        readterms-.->readtermsIter
        readtermsIter -.->define
        define-->fillunclass
        fillunclass-->readlith
        readlith-->mergelith
    end

    subgraph get_layers["Get (Hydro)geologic Layers"]
        direction TB

        subgraph uniquewells
            direction RL
            W0(("get_unique_wells()"))
            W1[(df)]

            W1-->W0
        end

        subgraph sortdf
            direction RL
            X0(("sort_dataframe()"))
            X1[(df)]
            X2{{sort_cols}}
            X3{{remove_nans}}

            X1-->X0
            X2-->X0
            X3-.->X0
        end

        subgraph alignrasters
            direction RL
            Y0(("align_rasters()"))
            Y1[(grids_unaligned)]
            Y2[(modelgrid)]

            Y1-->Y0
            Y2-->Y0
        end

        subgraph driftthick
            direction RL
            Z0(("get_drift_thick()"))
            Z1[(surface)]
            Z2[(bedrock)]
            Z3{{layers}}
            Z4{{plot}}

            Z1-->Z0
            Z2-->Z0
            Z3-.->Z0
            Z4-.->Z0
        end

        subgraph sampleraster
            direction RL
            AA0(("sample_raster_points()"))
            AA1[(raster)]
            AA2[(points_df)]
            AA3{{new_col}}

            AA1-->AA0
            AA2-->AA0
            AA3-.->AA0
        end

        sampleRasterIter[["Iterate sample_raster_points() as needed"]]

        subgraph layerdepths
            direction RL
            BB0(("get_layer_depths()"))
            BB1[(well_metadata)]
            BB2{{no_layers}}

            BB1-->BB0
            BB2-.->BB0
        end

        subgraph mergetables
            direction RL
            CC0(("merge_tables()"))
            CC1[(data_df)]
            CC2[(header_df)]
            CC3{{data_cols=None}}
            CC4{{header_cols=None}}
            CC5{{auto_pick_cols}}
            CC6{{"pd.merge(**kwargs)"}}

            CC1-->CC0
            CC2-->CC0
            CC3-.->CC0
            CC4-.->CC0
            CC5-.->CC0
            CC6-.->CC0
        end

        subgraph layertarg
            direction RL
            DD0(("layer_target_thick()"))
            DD1[(df)]
            DD2{{depth_top_col='TOP'}}
            DD3{{depth_bot_col='BOTTOM'}}

            DD1-->DD0
            DD2-.->DD0
            DD3-.->DD0
        end

        subgraph layerinterp
            direction RL
            EE0(("layer_interp()"))
            EE1[(points)]
            EE2[(grid)]
            EE3{{return_type='dataarray'}}
            EE4{{targetcol='TARG_THICK_PER'}}
            EE5{{"scipy.interpolate(**kwargs)"}}

            EE1-->EE0
            EE2-->EE0
            EE3-.->EE0
            EE4-.->EE0
            EE5-.->EE0  
        end

        uniquewells-->sortdf
        sortdf-->alignrasters
        alignrasters-->driftthick
        driftthick-->sampleraster
        sampleraster-->sampleRasterIter
        sampleRasterIter-->layerdepths
        layerdepths-->mergetables
        mergetables-->layertarg
        layertarg-->layerinterp

    end
    
    subgraph export["Export"]
        direction RL

        subgraph exportdf
            direction RL
            GG0(("export_dataframe()"))
            GG1[(df)]
            GG2{{out_dir}}
            GG3{{filename}}
            GG4{{datestamp=True}}

            GG1-->GG0
            GG2-->GG0
            GG3-->GG0
            GG4-.->GG0
        end

        subgraph exportdata
            direction RL
            HH0(("export_grids()"))
            HH1[(grid_data)]
            HH2{{out_path}}
            HH3{{file_id}}
            HH4{{filetype='tif'}}
            HH5{{variable_sep=True}}
            HH6{{date_stamp=True}}

            HH1-->HH0
            HH2-->HH0
            HH3-.->HH0
            HH4-.->HH0
            HH5-.->HH0
            HH6-.->HH0
        end

        exportdf~~~exportdata
    end

    setup-->read
    read-->readSA
    readSA-->Read_Grids
    Read_Grids-->clipdatatoSA
    clipdatatoSA-->clean
    clean-->classify
    classify-->get_layers
    get_layers-->export

```
