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
    
    # Next steps:
    # 1) For each profile/linestring, get segments
    # 2) For each segment, generate coordinates to use
    # 3) For each coordinate, get data for all layers
    # 4) Create a DataArray for each segment
    # 5) Plot all DataArrays in subplot
    # 6) Plot map of profile(s)
    
    return