    """The visualization module contains functions for creating plots and data vizualizations from outputs from w4h.run()"""
    
def plot_profile(dataset, profile=None, profile_direction=None, xcoord='x', ycoord='y', verbose=False, **kwargs):
    if profile is None:
        minX = np.min(dataset.coords[xcoord].values)
        midX = np.mean(dataset.coords[xcoord].values)
        maxX = np.max(dataset.coords[xcoord].values)

        minY = np.min(dataset.coords[ycoord].values)
        midY = np.mean(dataset.coords[ycoord].values)
        maxY = np.max(dataset.coords[ycoord].values)
        
        profileX = shapely.LineString([[minX, midY], [maxX, midY]])
        profileY = shapely.LineString([[midX, minY], [midX, maxY]])
        profile = [profileX, profileY]
        profile_direction = ['WE', 'SN']
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
    


    pLengths = []
    pSegLengths = []
    pSegOrients = []
    profileSamplePoints = []
    profileDicts = []
    # Iterate through each profile
    for p in profile:
        # Temporary lists for each segment/vertex
        segLengths = []
        segOrientations = []
        segSamplePoints = []

        # Iterate through each segment/vertex of the profile
        for i, vertex in enumerate(p.coords):
            if i==0:
                continue
            prevVertex = p.coords[i-1]

            segmentLineString = shapely.LineString([prevVertex, vertex])
            segXCoords = dataset.sel(x=slice(prevVertex[0], vertex[0])).coords['x'].values
            segYCoords = dataset.sel(x=slice(prevVertex[1], vertex[1])).coords['x'].values
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
        surfaceElevs = []
        bedrockElevs = []
        layerElevs = []
        modelData = []
        for point in segSamplePoints:
            xcoord, ycoord = point.xy
            surfaceElevs.append(dataset['Surface_Elevation'].sel(x=xcoord, y=ycoord, method='nearest').values.flatten)
            bedrockElevs.append(dataset['Bedrock_Elevation'].sel(x=xcoord, y=ycoord, method='nearest').values)
            modelData.append(dataset['Model_Layers'].sel(x=xcoord, y=ycoord, method='nearest').values)
            layerElevs.append(dataset.coords['layer_elevs'].sel(x=xcoord, y=ycoord, method='nearest').values)
        
        surfElevArr = np.array(surfaceElevs).flatten()
        bedElevArr = np.array(bedrockElevs).flatten()

        profileDicts.append({'Surface_Elevation': surfElevArr,
                            'Bedrock_Elevation': bedElevArr,
                            'Model_Layers': np.array(modelData).reshape((surfElevArr.shape[0], 9)),
                            'layer_elevs': np.array(layerElevs).reshape((surfElevArr.shape[0], 9))
                            })
            
    # Now just need to plot, got data

    return profileDicts