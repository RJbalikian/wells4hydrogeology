"""The visualization module contains functions for creating plots and data vizualizations from outputs from w4h.run()"""
import warnings

import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import pyproj
import shapely

def plot_cross_section(dataset, profile=None, profile_direction=None,
                       xcoord='x', ycoord='y',
                       mapped_variable='Depth_to_Bedrock',
                       cross_section_variable='Model_Layers',
                       surface_elevation_variable='Surface_Elevation',
                       bedrock_elevation_variable='Bedrock_Elevation',
                       layer_elevation_coordinate='layer_elevs',
                       show_layers=True, return_profile_dicts=False,
                       elev_unit='feet', convert_elevation_to=None,
                       title=None,
                       verbose=False, **kwargs):
    """Function to plot cross section profiles for datasets with properly configured coordinates and variables.
    This is intended to work "out of the box" with the xarray.Datasets output from w4h.run()

    Parameters
    ----------
    dataset : xarray.Dataset
        The xarray.Dataset with the proper data variables. Works "out of the box" with outputs from w4h.run().
    profile : None, shapely.Linestring, list of coordinates, or geopandas.GeoDataFrame, optional
        The profile(s) for which to create the cross sections.
        If None, by default creates one X profile and one Y profile in the middle of each dimension, by default None
    profile_direction : list of str, optional
        List of strings (list is same length as profile) indicating the direction to use for the profile. 
        If None, will be ['WE', 'SN'] to fit with profile=None defaults, by default None
    xcoord : str, optional
        Name of x coordinate, by default 'x'
    ycoord : str, optional
        Name of y coordinate, by default 'y'
    mapped_variable : str, optional
        Name of variable to show in main map, by default 'Depth_to_Bedrock'
    cross_section_variable : str, optional
        Name of variable to use for cross section profiles, by default 'Model_Layers'
    surface_elevation_variable : str, optional
        Variable to use for the surface elevation, by default 'Surface_Elevation'
    bedrock_elevation_variable : str, optional
        Variable to use for the bedrock elevation, by default 'Bedrock_Elevation'
    layer_elevation_coordinate : str, optional
        Coordinate name to use for the layer elevations. 
        This should be a non-indexed coordinate with the shape of the x, y, and layer coordinates, by default 'layer_elevs'
    show_layers : bool, optional
        Whether to plot the layer boundaries on the cross section, by default True
    return_profile_dicts : bool, optional
        Whether to return the profile dictionaries, rather than the matplotlib.Figure, by default False
    elev_unit : str, optional
        Unit of elevation for the elevation data, by default 'feet'
    convert_elevation_to : str, optional
        If None (default), does not convert elevation. Otherwise, will convert elevation to specified unit.
        Only conversion between 'ft' and 'meters' supported.
    title : str, optional
        Title to use for the output figure. If None, will be derived from variable names, by default None
    verbose : bool, optional
        Whether to print information about process to terminal, by default False

    Returns
    -------
    matplotlib.Figure
        Matplotlib.Figure instance is returned, unles return_profile_dicts is True.
        If return_profile_dicts=True, then a list of dicts with information about the profiles is returned.
    """
    
    
    if profile is None:
        firstX = dataset.coords[xcoord].values[0]
        midX = np.mean(dataset.coords[xcoord].values)
        lastX = dataset.coords[xcoord].values[-1]

        firstY = dataset.coords[ycoord].values[0]
        midY = np.mean(dataset.coords[ycoord].values)
        lastY = dataset.coords[ycoord].values[-1]
        
        profileX = shapely.LineString([[firstX, midY], [lastX, midY]])
        profileY = shapely.LineString([[midX, firstY], [midX, lastY]])
        profile = [profileX, profileY]
        #profile_direction = ['WE', 'SN']
    elif isinstance(profile, gpd.GeoDataFrame):
        profile = profile.geometry.to_list()
    elif isinstance(profile, shapely.Geometry):
        profile = [profile]
    elif isinstance(profile, (list, tuple)):
        if not all(isinstance(item, (shapely.Geometry, gpd.GeoDataFrame)) for item in profile):
            print("If using a list or tuple for profile parameter, all items must be a shapely Geometry or geopandas GeoDataFrame")
            return
        newProfiles = profile.copy()
        for p in profile:
            if isinstance(p, gpd.GeoDataFrame):
                profileExt = p.geometry.to_list()
                newProfiles.extend(profileExt)
        profile = newProfiles
    else:
        print(f"\tThe profile parameter must shapely Geometry or geopandas GeoDataFrame or a list of those objects")
        return 
    
    if profile_direction is None:
        profile_direction = []
        for p in profile:
            firstCoord = p.coords[0]
            lastCoord = p.coords[-1]
            xRange = abs(lastCoord[0] - firstCoord[0])
            yRange = abs(lastCoord[1] - firstCoord[1])
            
            if yRange > xRange:
                if lastCoord[1] > firstCoord[1]:
                    profile_direction.append("SN")
                else:
                    profile_direction.append("NS")
            else:
                if lastCoord[0] > firstCoord[0]:
                    profile_direction.append("WE")
                else:
                    profile_direction.append("EW")

    if len(profile_direction) != len(profile):
        print(f"\tprofile and profile_direction must be the same length, but len(profile)={len(profile)} and len(profile_direction)={len(profile_direction)}")
        return


    if convert_elevation_to is not None:
        if elev_unit.lower() in ['feet', 'foot', 'ft', 'f']:
            if str(convert_elevation_to).lower() in ['meters', 'metres', 'meter', 'metre', 'mtr', 'm']:
                dataset[mapped_variable] = dataset[mapped_variable] * 0.3048
                dataset[surface_elevation_variable] = dataset[surface_elevation_variable] * 0.3048
                dataset[bedrock_elevation_variable] = dataset[bedrock_elevation_variable] * 0.3048
                dataset[layer_elevation_coordinate] = dataset[layer_elevation_coordinate] * 0.3048
                elev_unit = 'meters'
        elif elev_unit.lower() in ['meters', 'metres', 'meter', 'metre', 'mtr', 'm']:
            if str(convert_elevation_to).lower() in ['feet', 'foot', 'ft', 'f']:
                dataset[mapped_variable] = dataset[mapped_variable] * 3.2808399
                dataset[surface_elevation_variable] = dataset[surface_elevation_variable] * 3.2808399
                dataset[bedrock_elevation_variable] = dataset[bedrock_elevation_variable] * 3.2808399
                dataset[layer_elevation_coordinate] = dataset[layer_elevation_coordinate] * 3.2808399
                elev_unit = 'feet'

    profileDicts = []
    # Iterate through each profile
    for p_i, p in enumerate(profile):
        # Temporary lists for each segment/vertex
        segLengths = []
        segOrientations = []
        segSamplePoints = []
        surfaceElevs = []
        bedrockElevs = []
        layerElevs = []
        modelData = []
        
        # Iterate through each segment/vertex of the profile
        for i, vertex in enumerate(p.coords):
            if i == 0:
                continue
            prevVertex = p.coords[i-1]

            segmentLineString = shapely.LineString([prevVertex, vertex])
            segXCoords = dataset.sel(x=slice(prevVertex[0], vertex[0])).coords['x'].values
            segYCoords = dataset.sel(x=slice(prevVertex[1], vertex[1])).coords['y'].values

            segXCoordsDists = np.diff(segXCoords)
            segYCoordsDists = np.diff(segYCoords)
            segLengths.append(segmentLineString.length)

            # Get distance of each segment
            y1 = vertex[1]
            y0 = prevVertex[1]
            x1 = vertex[0]
            x0 = prevVertex[0]
            yDist = y1-y0
            xDist = x1-x0

            # Get orientation of segment
            segOrient='WE'
            if abs(yDist) > abs(xDist):
                segOrient = 'SN'
                if yDist < 0:
                    segOrient = 'NS'
            else:
                if xDist < 0:
                    segOrient = 'EW'
            segOrientations.append(segOrient)
            yOrients = ['NS', "SN", 'NORTH', 'SOUTH', 'N', 'S', 'LATITUDE', 'NORTHING', "SOUTHING"]
            xOrients = ['WE', "EW", 'EAST', 'WEST', 'E', 'W', "LONGITUDE", "EASTING", "WESTING"]

            # Get coordinates to sample for profile
            if segOrient.upper() in xOrients:
                vertexIndex = 0
                otherVertextInd = 1
                segCoords = segXCoords
                segCoordsDists = segXCoordsDists 

            else: # if oriented more NS
                vertexIndex = 1
                otherVertextInd = 0
                segCoords = segYCoords
                segCoordsDists = segYCoordsDists 

            # For each segment, get each xy value as shapely point at each existing coordinate value
            subsegSamplePoints = []
            for i, segmentSegment in enumerate(segCoordsDists):
                currSubsegDistance = np.sum(segCoordsDists[:i+1])
                subsegSamplePoints.append(segmentLineString.line_interpolate_point(distance=currSubsegDistance))
            segSamplePoints.extend(subsegSamplePoints)

        profileLineString = shapely.LineString(segSamplePoints)
        for point in segSamplePoints:
            xcoord, ycoord = point.xy
            surfaceElevs.append(dataset[surface_elevation_variable].sel(x=xcoord, y=ycoord, method='nearest').values)
            bedrockElevs.append(dataset[bedrock_elevation_variable].sel(x=xcoord, y=ycoord, method='nearest').values)
            modelData.append(dataset[cross_section_variable].sel(x=xcoord, y=ycoord, method='nearest').values)
            layerElevs.append(dataset.coords[layer_elevation_coordinate].sel(x=xcoord, y=ycoord, method='nearest').values)
    
        surfElevArr = np.array(surfaceElevs).flatten()
        bedElevArr = np.array(bedrockElevs).flatten()

        profileDicts.append({'XY': shapely.LineString(segSamplePoints),
                             'Surface_Elevation': surfElevArr,
                             'Bedrock_Elevation': bedElevArr,
                             'xSection_Data': np.array(modelData).reshape((surfElevArr.shape[0], 9)),
                             'layer_elevs': np.array(layerElevs).reshape((surfElevArr.shape[0], 9)),
                             'profile_direction': profile_direction[p_i],
                             'elev_unit': elev_unit
                             })

    # Now just need to plot, already got data
    # Set up entire figure
    noSubplots = len(profileDicts)
    xSecSPList = []
    for letter in range(ord('A'), ord('A')+noSubplots):
        i = letter - ord('A')
        subplotName = 'XSEC'+chr(letter)
        profileDicts[i]['ProfileName'] = chr(letter)
        profileDicts[i]['SubplotName'] = subplotName

        xSecSPList.append([subplotName])

    subplotMosaicLists = [["MAP"],
                          ["MAP"],
                          ["MAP"]]
    subplotMosaicLists.extend(xSecSPList)

    # Get Figure/Axes sizing
    figH = 11
    figW = 17
    if len(subplotMosaicLists)<2:
        pass
    elif len(subplotMosaicLists) < 4:
        figH = 17
        figW = 11
    else:
        figW = 11
        figH = 11 + len(subplotMosaicLists) * 1.5
    plt.rcParams['figure.figsize'] = (figW, figH)
    fig, ax = plt.subplot_mosaic(subplotMosaicLists, layout='tight')

    if title is None:
        mVarTitle = mapped_variable.replace('_', ' ').title()
        xSecTitle = cross_section_variable.replace('_', ' ').title()
        titleStr = f'Map of {mVarTitle} with profile(s) showing {xSecTitle}'
    else:
        titleStr = title

    # Code to make map
    colorbar_kwargs = {'label': f"{mVarTitle} [{elev_unit}]"}
    dataset[mapped_variable].plot(ax=ax['MAP'],
                                cbar_kwargs=colorbar_kwargs)
    
    ax['MAP'].ticklabel_format(style='plain')
    ax['MAP'].set_aspect('equal')
    ax['MAP'].set_title(titleStr)

    proj = pyproj.CRS.from_user_input(dataset.spatial_ref.attrs['crs_wkt'])
    xUnit = 'Deg.'
    yUnit = 'Deg.'
    xName = 'Longitude'
    yName = 'Latitude'

    for a in proj.axis_info:
        if hasattr(a, 'direction'):
            if a.direction.upper() in xOrients:
                xUnit = str(a.unit_name)
                xName = str(a.name)
            else:
                yUnit = str(a.unit_name)
                yName = str(a.direction)

    xCoordLabel = f"{xName.upper()} [{xUnit.upper()}]"
    yCoordLabel = f"{yName.upper()} [{yUnit.upper()}]"
    ax['MAP'].set_xlabel(xCoordLabel)
    ax['MAP'].set_ylabel(yCoordLabel)

    # Code for plotting each profile
    for profile in profileDicts:
        # Set up plot
        currSubP = profile['SubplotName']
        plt.sca(ax[currSubP])

        plt.rcParams['xtick.top'] = True
        plt.rcParams["xtick.labeltop"] = True
        plt.rcParams["ytick.right"] = True
        plt.rcParams["ytick.labelright"] = True
        ax[currSubP].ticklabel_format(style='plain')

        xArray, yArray = profile["XY"].xy

        xArray = np.array(xArray.tolist())
        yArray = np.array(yArray.tolist())
        layer_elevs = profile['layer_elevs']
        xSection_Data = profile["xSection_Data"]
        br_elev = profile['Bedrock_Elevation']
        surf_elev = profile['Surface_Elevation']
        # Find rows in layer_elevs that contain any NaN values
        valid_rows = ~np.isnan(layer_elevs).any(axis=1)

        # Filter xArray, layer_elevs, and xSection_Data based on valid rows
        xArray_filtered = xArray[valid_rows]
        yArray_filtered = yArray[valid_rows]
        br_elevs_filtered = br_elev[valid_rows]
        surf_elevs_filtered = surf_elev[valid_rows]
        layer_elevs_filtered = layer_elevs[valid_rows, :]
        xSection_Data_filtered = xSection_Data[valid_rows, :]


        # Plot profile on map
        ax['MAP'].plot(xArray_filtered, yArray_filtered,
                       c='k', linewidth=2,#linestyle='dashed',
                       path_effects=[pe.withStroke(linewidth=3.5,
                                                   foreground="w")])
        
        # Add profile annotation
        ax['MAP'].text(xArray_filtered[0], yArray_filtered[0],
                       profile['ProfileName'],
                       path_effects=[pe.withStroke(linewidth=4, foreground="w")])
        
        ax['MAP'].text(xArray_filtered[-1], yArray_filtered[-1],
                       profile['ProfileName']+"'",
                       path_effects=[pe.withStroke(linewidth=4, foreground="w")])
        

        if profile['profile_direction'].upper() in xOrients:
            X = np.tile(xArray_filtered, (layer_elevs_filtered.shape[1], 1)).T
            coords = xArray_filtered
            coordLabel = xCoordLabel
        else:
            X = np.tile(yArray_filtered, (layer_elevs_filtered.shape[1], 1)).T
            coords = yArray_filtered
            coordLabel = yCoordLabel

        Y = layer_elevs_filtered

        # Plot data
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pcm = ax[currSubP].pcolormesh(X, Y, xSection_Data_filtered, cmap='Oranges', alpha=0.8)
            #fig.colorbar(pcm, ax=ax['MAP'])
        except Exception as e:
            print(f'colormesh didnt work: {e}')

        # Show model layers?
        if show_layers:
            ax[currSubP].plot(coords, layer_elevs_filtered,
                              c='k', linewidth=0.5)
            
        # Show cell divisions (vertical)?
        if len(coords) < 1000:
            ax[currSubP].vlines(coords,
                                ymin=br_elevs_filtered,
                                ymax=surf_elevs_filtered,
                                linewidths=0.25, colors='k')
            
        
        minBR = np.nanmin(br_elevs_filtered)
        maxSurf = np.nanmax(surf_elevs_filtered)
        
        stopVal = np.nanmax(np.subtract(surf_elev, br_elev))/(9/1.25)
        stepVal = stopVal//10
        ax[currSubP].fill_between(coords, br_elevs_filtered, br_elevs_filtered-(stopVal), facecolor='white', linewidth=2)
        for i in np.arange(0, stopVal, stepVal):
            ax[currSubP].fill_between(coords, br_elevs_filtered, br_elevs_filtered-(stopVal-i), facecolor='purple', linewidth=2, alpha=0.1)

        ax[currSubP].plot(coords, br_elevs_filtered, c='purple', linewidth=2)
        ax[currSubP].plot(coords, surf_elevs_filtered, c='k', linewidth=2)
        
        # Set y limits
        yRange = abs(maxSurf - (minBR-stopVal))
        yHi = maxSurf + (yRange * 0.15) # Need more space on top for annotation
        yLo = (minBR-stopVal) - (yRange * 0.05)
        ax[currSubP].set_ylim([yLo, yHi])

        ax[currSubP].annotate(text=profile['ProfileName'],
                              xy=(0.01, 0.98),
                              xycoords='axes fraction',
                              verticalalignment='top',
                              horizontalalignment='left',
                              path_effects=[pe.withStroke(linewidth=4,
                                                          foreground="w")])
        
        ax[currSubP].annotate(text=profile['ProfileName']+"'",
                              xy=(0.99, 0.98),
                              xycoords='axes fraction',
                              verticalalignment='top',
                              horizontalalignment='right',
                              path_effects=[pe.withStroke(linewidth=4,
                                                          foreground="w")])

        ax[currSubP].set_xlabel(coordLabel)
        ax[currSubP].set_ylabel(f"Elevation [{profile['elev_unit']}]")
        
    plt.show()

    if return_profile_dicts:
        return profileDicts

    return fig