"""The visualization module contains functions for creating plots and data vizualizations from outputs from w4h.run()"""
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import shapely

def plot_cross_section(dataset, profile=None, profile_direction=None, xcoord='x', ycoord='y', verbose=False, **kwargs):
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
            if i==0:
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
            yOrients = ['NS', "SN"]
            xOrients = ['WE', "EW"]

            # Get coordinates to sample for profile
            if segOrient in xOrients:
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
            surfaceElevs.append(dataset['Surface_Elevation'].sel(x=xcoord, y=ycoord, method='nearest').values)
            bedrockElevs.append(dataset['Bedrock_Elevation'].sel(x=xcoord, y=ycoord, method='nearest').values)
            modelData.append(dataset['Model_Layers'].sel(x=xcoord, y=ycoord, method='nearest').values)
            layerElevs.append(dataset.coords['layer_elevs'].sel(x=xcoord, y=ycoord, method='nearest').values)
        
        surfElevArr = np.array(surfaceElevs).flatten()
        bedElevArr = np.array(bedrockElevs).flatten()

        profileDicts.append({'XY': shapely.LineString(segSamplePoints),
                            'Surface_Elevation': surfElevArr,
                            'Bedrock_Elevation': bedElevArr,
                            'Model_Layers': np.array(modelData).reshape((surfElevArr.shape[0], 9)),
                            'layer_elevs': np.array(layerElevs).reshape((surfElevArr.shape[0], 9)),
                            'profile_direction': profile_direction[p_i]
                            })
            
    # Now just need to plot, got data
    # FROM HERE ON DOWN DOES NOT WORK PERFECTLY, NEED TO ADJUST FOR DIFFERENT ORIENTATIONS
    show_layers=True
    for profile in profileDicts:
        xArray, yArray = profile["XY"].xy
        xArray = np.array(xArray.tolist())
        yArray = np.array(yArray.tolist())
        layer_elevs = profile['layer_elevs']
        Model_Layers = profile["Model_Layers"]
        br_elev = profile['Bedrock_Elevation']
        surf_elev = profile['Surface_Elevation']
        # Find rows in layer_elevs that contain any NaN values
        valid_rows = ~np.isnan(layer_elevs).any(axis=1)

        # Filter xArray, layer_elevs, and Model_Layers based on valid rows
        xArray_filtered = xArray[valid_rows]
        layer_elevs_filtered = layer_elevs[valid_rows, :]
        Model_Layers_filtered = Model_Layers[valid_rows, :]

        X = np.tile(xArray_filtered, (layer_elevs_filtered.shape[1], 1)).T
        Y = layer_elevs_filtered

        # Set up plot
        plt.rcParams['figure.figsize'] = (30, 8)
        plt.rcParams['xtick.top'] = True
        plt.rcParams["xtick.labeltop"] = True
        plt.rcParams["ytick.right"] = True
        plt.rcParams["ytick.labelright"] = True


        # Plot data
        try:
            pcm = plt.pcolormesh(X, Y, Model_Layers_filtered, cmap='Oranges', alpha=0.8)
        except Exception as e:
            print('colormesh didnt work')

        if show_layers:
            plt.plot(xArray, profile['layer_elevs'], c='k', linewidth=0.5)
        if len(xArray)< 1000:
            plt.vlines(xArray, ymin=profile['Bedrock_Elevation'], ymax=profile['Surface_Elevation'], linewidths=0.25, colors='k')
        minBR = np.nanmin(profile['Bedrock_Elevation'])
        stopVal=np.nanmax(np.subtract(surf_elev, br_elev))/(9/1.25)
        stepVal = stopVal//10
        plt.fill_between(xArray, br_elev, br_elev-(stopVal), facecolor='white', linewidth=3)
        for i in np.arange(0, stopVal, stepVal):
            plt.fill_between(xArray, br_elev, br_elev-(stopVal-i), facecolor='purple', linewidth=3, alpha=0.1)

        plt.plot(xArray, br_elev, c='purple', linewidth=3)
        plt.plot(xArray, surf_elev, c='k', linewidth=3)
        plt.show()

    return profileDicts