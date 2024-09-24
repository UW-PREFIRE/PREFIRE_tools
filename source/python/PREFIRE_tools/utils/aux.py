import numpy as np
import netCDF4 as n
import pyproj
import pygeos
import datetime

from PREFIRE_PRD_GEN.file_read import get_PREFIRE_Lx_field


def construct_tirs_geo(tirs_fpath, artp, polygons=True):
    """
    From the input PREFIRE product file, read the time of each TIRS
    field-of-view (FOV) and the shape of the dimensions specifying the number
    of FOVs.
    
    Optionally, create polygons representing the outline of each TIRS
    FOV.

    Parameters
    ----------
    tirs_fpath : str
        Path to a PREFIRE product file containing a Geometry group.
    artp : 3-tuple (str, int, int)
        (optional) A 3-tuple containing the dimension name to subset, and start
         and stop indices (NumPy indexing convention) in the given granule
    polygons : bool, optional
        Whether to create a polygon representing the outline of each TIRS FOV.
        The default is True.

    Raises
    ------
    Exception
        Error: TIRS FOV split into more than 2 polygons.

    Returns
    -------
    tirs_times : list
        Time of each TIRS FOV (type: datetime.datetime)
    tirs_shape : tuple
        Shape of the dimensions specifying the number of FOVs (atrack, xtrack)
    tirs_lons : numpy.ndarray
        Longitude of each FOV's centroid [deg_E]
    tirs_lats : numpy.ndarray
        Latitude of each FOV's centroid [deg_N]
    tirs_polylist : list, optional
        Polygon outline(s) for each TIRS FOV. Each list item is a list of
        pygeos polygons that has length 2 if the TIRS polygon crosses the
        antimeridian.
    """
    
    with n.Dataset(tirs_fpath) as nc:
        tirs_lons = get_PREFIRE_Lx_field(nc, "Geometry", "longitude", artp)
        tirs_shape = tirs_lons.shape
        tirs_lats = get_PREFIRE_Lx_field(nc, "Geometry", "latitude", artp)

        # Reshape TIRS FOV vertex lats/lons from 3D to 2D arrays
        # - 3D dims: (atrack, xtrack, FOV_vertices)
        # - 2D dims: (atrack * xtrack, FOV_vertices)
        tmpB = get_PREFIRE_Lx_field(nc, "Geometry", "vertex_longitude", artp)
        vlons_reshp = np.reshape(tmpB, (np.size(tirs_lons), 4))
        tmpB = get_PREFIRE_Lx_field(nc, "Geometry", "vertex_latitude", artp)
        vlats_reshp = np.reshape(tmpB, (np.size(tirs_lons), 4))

        tmpC = get_PREFIRE_Lx_field(nc, "Geometry", "time_UTC_values", artp)
        times_split = np.array(tmpC)  #? Is this np.array() call needed?

    # TIRS times are only provided for each along-track swath. Duplicate
    # times for each cross-track scene to match length of TIRS polygon list.
    tirs_times = []
    for ai in range(tirs_shape[0]):
        t_array = times_split[ai]
        for xi in range(tirs_shape[1]):
            tirs_times.append(TIRS_L1B_time_array_to_dt(t_array))
    
    if polygons:
        tirs_polylist = []
        for vlons_fov, vlats_fov in zip(vlons_reshp, vlats_reshp):
            vlons_closed = np.append(vlons_fov, vlons_fov[0])
            vlats_closed = np.append(vlats_fov, vlats_fov[0])
            
            lon_diffs = [lon - vlons_closed for lon in vlons_closed]
            # Crosses antimeridian (TIRS scenes should never cross a pole)
            if np.max(np.abs(lon_diffs)) > 180:
                split_polys = []
                poly_pts_unsplit = generate_gc_poly_bounds(
                    vlons_closed, vlats_closed
                    )                
                split_polys_pts = poly_antimeridian_split(poly_pts_unsplit)
                for poly_pts in split_polys_pts:
                    split_polys.append(pygeos.polygons(poly_pts))
                
                # Raise exception if TIRS scene split into more than 2 polygons.
                # This would indicate a bug in poly_antimeridian_split.
                if len(split_polys) > 2:
                    raise Exception(
                        'Error: TIRS scene split into more than 2 polygons.'
                        )
                    
                tirs_polylist.append(split_polys)
            else:
                poly = pygeos.polygons(
                    [[lon,lat] for lon,lat in zip(vlons_closed, vlats_closed)]
                    )                
                tirs_polylist.append([poly])
            
        return (tirs_times, tirs_shape, tirs_lons, tirs_lats, tirs_polylist)
    
    else:
        return (tirs_times, tirs_shape, tirs_lons, tirs_lats)


def generate_gc_poly_bounds(vert_lons_ccw, vert_lats_ccw, npts=200):   
    """
    Create a list of polygon edge points that are evenly-spaced along great 
    circle lines between the polygon vertices.
    
    Parameters
    ----------
    vert_lons_ccw : numpy.ndarray
        Longitudes of polygon vertices in counter-clockwise order.
    vert_lats_ccw : numpy.ndarray
        Latitudes of polygon vertices in counter-clockwise order.
    npts : int, optional
        Number of equally spaced points along the geodesic between successive
        vertices. The default is 200.
    
    Returns
    -------
    poly_gc_pts : list
        List of points defining the polygon boundaries.
    
    """
    geod = pyproj.Geod(ellps='WGS84')
    
    poly_gc_pts = []
    for vert_ix in range(len(vert_lons_ccw)-1):
        poly_gc_pts.extend(geod.npts(vert_lons_ccw[vert_ix], 
                                     vert_lats_ccw[vert_ix], 
                                     vert_lons_ccw[vert_ix+1],
                                     vert_lats_ccw[vert_ix+1],
                                     npts,
                                     initial_idx=0,
                                     terminus_idx=1))
        # geod.npts doesn't include end points, so append manually to
        # ensure polygon edges aren't "rounded"
        poly_gc_pts.append((vert_lons_ccw[vert_ix+1], vert_lats_ccw[vert_ix+1]))
    
    return poly_gc_pts


def check_pole_antimeridian_crosses(poly_bound_pts):
    """
    "Walk" along the outer boundary points of a polygon.

    Polygons that cross neither the antimeridian nor a pole will have no sudden
    discontinuities in longitude values.
    
    Polygons that cross the antimeridian but not a pole will have two discontinuities.
    
    Polygons that cross a single pole will have one discontinuity. (JPSS polygons will never cross 
    both poles; TIRS polygons will never cross a pole at all.)
    
    This function works with both TIRS and JPSS polygons.
    
    Parameters
    ----------
    poly_bound_pts : list
        Coordinates of polygon boundary points, in counter-clockwise order (i.e. the interior
        of the polygon is on your left-hand side as you are walking the boundary).
    
    Returns
    -------
    num_crosses : int
        Integer with a value of 0 (no crosses), 1 (pole cross only), or 2 (2 antimeridian 
        crosses, therefore no pole cross).
    
    """
        
    lon_diffs = [poly_bound_pts[0][0] - poly_bound_pts[-1][0]]
    for i in range(len(poly_bound_pts)-1):
        lon_diff = poly_bound_pts[i+1][0] - poly_bound_pts[i][0]
        lon_diffs.append(np.abs(lon_diff))
    
    discts = np.where(np.array(lon_diffs) > 180)
    num_crosses = np.shape(discts)[1]
    
    return num_crosses


def poly_antimeridian_split(poly_bound_pts):
    """
    Split polygons that cross the antimeridian into two polygons that have the
    -180 longitude line as one edge. An exception is raised if the longitude
    and latitude vertices passed to the function define a polygon that either
    (a) doesn't cross the antimeridian, or (b) crosses a pole.
    
    This function works with both TIRS and JPSS polygons in the "simple" case
    of a polygon with no convex edges that cross the antimeridian more than
    once. The "convex" edge case only occurs for JPSS polygons and is handled
    by the _convex_split helper function.
    
    ** NOTE: To create "precise" polygons split across the antimeridian, the
    boundary points of the polygon should be closely spaced, like the TIRS and
    VIIRS polygons created by construct_tirs_geo and
    create_VIIRS_poly_from_outer_pts, respectively. If only the four corner
    points of a TIRS polygon are used, the split polygon shape will be
    imprecise because the mean latitude of the widely-spaced corner points will
    be used as the latitude of the "split" point along the antimeridian.
    
    Parameters
    ----------
    poly_bound_pts : list
        Coordinates of polygon boundary points, in counter-clockwise order
        (i.e. the interior of the polygon is on your left-hand side as you are
         walking the boundary).
    
    Returns
    -------
    split_polys_list : list
        Each list entry is a list of points defining the boundaries of an
        individual polygon produced by the split.
    
    """

    # Verify that the polygon crosses the date line but not a pole.
    # Raise exception if num_crosses is 0 (no crossing) or 1 (pole crossing).
    num_crosses = check_pole_antimeridian_crosses(poly_bound_pts)
    if num_crosses == 0:
        raise Exception("Polygon does not cross antimeridian.")
    elif num_crosses == 1:
        raise Exception("Polygon crosses a pole. Can't split along "
                        "antimeridian.")
    
    # Walk through all points along the polygon edge to find where polygon will
    # be split. When a discontinuity is found, mark the list index and the type 
    # of discontinuity.
    # - split type 1: walking from eastern hemisphere into western hemisphere
    # - split type 2: walking from western hemisphere into eastern hemisphere

    split_ixs = []    
    split_types = []    
    
    # First, handle the case where antimeridian crossing happens between the 
    # last and first point in the list
    if poly_bound_pts[0][0] - poly_bound_pts[-1][0] < -180:
        split_ixs.append(-1)
        split_types.append(1)
    elif poly_bound_pts[0][0] - poly_bound_pts[-1][0] > 180:
        split_ixs.append(-1)
        split_types.append(2)
    else:
        pass
    
    # All other antimeridian crossing cases (occurring at any other point in
    # list)
    for i in range(len(poly_bound_pts)-1):
        lon_diff = poly_bound_pts[i+1][0] - poly_bound_pts[i][0]

        if lon_diff < -180:
            split_ixs.append(i)
            split_types.append(1)
        elif lon_diff > 180:
            split_ixs.append(i)
            split_types.append(2)
        else:
            pass
    
    # "Simple" split case: polygon does not have any convex edges that cross
    # the antimeridian more than once
    if num_crosses == 2:
        split_polys_list = _simple_split(
            poly_bound_pts, split_ixs, split_types)
        
    # "Convex" split case: polygon has a convex edge that crosses the
    # antimeridian more than once
    elif num_crosses == 4:
        split_polys_list = _convex_split(
            poly_bound_pts, split_ixs, split_types)
    
    return split_polys_list

def _simple_split(poly_bound_pts, split_ixs, split_types):
    """
    Helper to poly_antimeridian_split that handles the case where the polygon
    does not have any convex edges that cross the antimeridian more than once.
    
    Two polygons are created in the "simple" split case, and one edge of each
    polygon is created by drawing a boundary directly along the -180 / 180
    longitude line.

    Parameters
    ----------
    poly_bound_pts : list
        Coordinates of polygon boundary points, in counter-clockwise order
        (i.e. the interior of the polygon is on your left-hand side as you are
         walking the boundary).
    split_ixs : list
        List indices into poly_bound_pts where the polygon crosses the
        antimeridian.
    split_types : list
        Split type (1=EH to WH; 2=WH to EH) at each split_ix.

    Returns
    -------
    split_polys_list : list
        Each list entry is a list of points defining the boundaries of an
        individual polygon produced by the split.

    """
    
    # The latitudes of these edge points are calculated as the mean latitude of
    # the two points on either side of the split.
    split_lat_1 = (poly_bound_pts[split_ixs[0]][1] + \
                   poly_bound_pts[split_ixs[0]+1][1]) / 2
    split_lat_2 = (poly_bound_pts[split_ixs[1]][1] + \
                   poly_bound_pts[split_ixs[1]+1][1]) / 2

    # First, handle the case where antimeridian crossing happens between the 
    # last and first point in the list
    if split_ixs[0] == -1:
        if split_types[0] == 1:
            poly_1_pts = [(-180,split_lat_1)]
            poly_1_pts.extend(poly_bound_pts[:split_ixs[1]+1])
            poly_1_pts.append((-180,split_lat_2))
            poly_2_pts = [(180,split_lat_2)]
            poly_2_pts.extend(poly_bound_pts[split_ixs[1]+1:])
            poly_2_pts.append((180,split_lat_1))
        elif split_types[0] == 2:
            poly_1_pts = [(180,split_lat_1)]
            poly_1_pts.extend(poly_bound_pts[:split_ixs[1]+1])
            poly_1_pts.append((180,split_lat_2))
            poly_2_pts = [(-180,split_lat_2)]
            poly_2_pts.extend(poly_bound_pts[split_ixs[1]+1:])
            poly_2_pts.append((-180,split_lat_1))
            
    # All other antimeridian crossings (occurring at any other point in list)
    else:
        poly_1_pts = poly_bound_pts[:split_ixs[0]+1]
        if split_types[0] == 1:
            poly_1_pts.append((180,split_lat_1))
            poly_2_pts = [(-180,split_lat_1)]
            poly_2_pts.extend(poly_bound_pts[split_ixs[0]+1:split_ixs[1]+1])
            poly_2_pts.append((-180,split_lat_2))
            poly_1_pts.append((180,split_lat_2))
        elif split_types[0] == 2:
            poly_1_pts.append((-180,split_lat_1))
            poly_2_pts = [(180,split_lat_1)]
            poly_2_pts.extend(poly_bound_pts[split_ixs[0]+1:split_ixs[1]+1])
            poly_2_pts.append((180,split_lat_2))
            poly_1_pts.append((-180,split_lat_2))
        
        poly_1_pts.extend(poly_bound_pts[split_ixs[1]+1:])
    
    split_polys_list = [poly_1_pts, poly_2_pts]
    
    return split_polys_list
    
def _convex_split(poly_bound_pts, split_ixs, split_types):
    """
    Helper to poly_antimeridian_split that handles the case where the polygon
    has a convex edge that crosses the antimeridian more than once.
    
    Three polygons are created in the "convex" split case:
    - Two triangles with one edge formed by the antimeridian
    - One six-sided polygon with two edges formed by the antimeridian

    Parameters
    ----------
    poly_bound_pts : list
        Coordinates of polygon boundary points, in counter-clockwise order
        (i.e. the interior of the polygon is on your left-hand side as you are
         walking the boundary).
    split_ixs : list
        List indices into poly_bound_pts where the polygon crosses the
        antimeridian.
    split_types : list
        Split type (1=EH to WH; 2=WH to EH) at each split_ix.

    Returns
    -------
    split_polys_list : list
        Each list entry is a list of points defining the boundaries of an
        individual polygon produced by the split.

    """

    # The latitudes of these edge points are calculated as the mean latitude of
    # the two points on either side of the split.
    split_lat_1 = (poly_bound_pts[split_ixs[0]][1] + \
                   poly_bound_pts[split_ixs[0]+1][1]) / 2
    split_lat_2 = (poly_bound_pts[split_ixs[1]][1] + \
                   poly_bound_pts[split_ixs[1]+1][1]) / 2
    split_lat_3 = (poly_bound_pts[split_ixs[2]][1] + \
                   poly_bound_pts[split_ixs[2]+1][1]) / 2
    split_lat_4 = (poly_bound_pts[split_ixs[3]][1] + \
                   poly_bound_pts[split_ixs[3]+1][1]) / 2
        
    # Four possible cases of split types: "bowl" formed by the concave edge and
    # the antimeridian is either located in the EH or the WH (with two
    # triangles in the opposite hemisphere from the bowl), and orbit is either
    # ascending or descending
    # - This assumes that the polygon was created with points oriented
    #   counterclockwise, *starting at the top left corner* of the polygon
    # - The bottom left corner of the polygon is defined relative to the
    #   satellite's orbit, so for descending orbits, the orbit-relative
    #   top left corner is the bottom right corner of the polygon when viewed
    #   on a map
        
    # First, handle the case where antimeridian crossing happens between the 
    # last and first point in the list
    if split_ixs[0] == -1:
        # Case 1: EH bowl / WH triangles, descending orbit
        if split_lat_1 < split_lat_2 < split_lat_3 < split_lat_4:
            poly_1_pts = [(-180,split_lat_1)]
            poly_1_pts.extend(poly_bound_pts[:split_ixs[1]+1])
            poly_1_pts.append((-180,split_lat_2))
            poly_2_pts = [(180,split_lat_2)]
            poly_2_pts.extend(poly_bound_pts[split_ixs[1]+1:split_ixs[2]+1])
            poly_2_pts.append((180,split_lat_3))
            poly_3_pts = [(-180,split_lat_3)]
            poly_3_pts.extend(poly_bound_pts[split_ixs[2]+1:split_ixs[3]+1])
            poly_3_pts.append((-180,split_lat_4))
            poly_2_pts.append((180,split_lat_4))
            poly_2_pts.extend(poly_bound_pts[split_ixs[3]+1:split_ixs[0]+1])
            poly_2_pts.append((180,split_lat_1))
        
        # Case 2: WH bowl / EH triangles, ascending orbit
        elif split_lat_4 < split_lat_3 < split_lat_2 < split_lat_1:
            poly_1_pts = [(180,split_lat_1)]
            poly_1_pts.extend(poly_bound_pts[:split_ixs[1]+1])
            poly_1_pts.append((180,split_lat_2))
            poly_2_pts = [(-180,split_lat_2)]
            poly_2_pts.extend(poly_bound_pts[split_ixs[1]+1:split_ixs[2]+1])
            poly_2_pts.append((-180,split_lat_3))
            poly_3_pts = [(180,split_lat_3)]
            poly_3_pts.extend(poly_bound_pts[split_ixs[2]+1:split_ixs[3]+1])
            poly_3_pts.append((180,split_lat_4))
            poly_2_pts = [(-180,split_lat_4)]
            poly_2_pts.extend(poly_bound_pts[split_ixs[3]+1:split_ixs[0]+1])
            poly_2_pts.append((-180,split_lat_1))
            
        # (Don't think either case 3 [EH bowl / WH triangles, ascending orbit]
        # or case 4 [WH bowl / EH triangles, descending orbit] is possible,
        # assuming polygons are built by starting at the top left corner)          
    
    else:
        poly_1_pts = poly_bound_pts[:split_ixs[0]+1]
        
        # Case 1: EH bowl / WH triangles, ascending orbit
        if split_lat_1 < split_lat_2 < split_lat_3 < split_lat_4:
            poly_1_pts.append((180,split_lat_1))
            poly_2_pts = [(-180,split_lat_1)]
            poly_2_pts.extend(poly_bound_pts[split_ixs[0]+1:split_ixs[1]+1])
            poly_2_pts.append((-180,split_lat_2))
            poly_1_pts.append((180,split_lat_2))
            poly_1_pts.extend(poly_bound_pts[split_ixs[1]+1:split_ixs[2]+1])
            poly_1_pts.append((180,split_lat_3))
            poly_3_pts = [(-180,split_lat_3)]
            poly_3_pts.extend(poly_bound_pts[split_ixs[2]+1:split_ixs[3]+1])
            poly_3_pts.append((-180,split_lat_4))
            poly_1_pts.append((180,split_lat_4))
            
        # Case 2: WH bowl / EH triangles, ascending orbit
        elif split_lat_3 < split_lat_2 < split_lat_1 < split_lat_4:
            poly_1_pts.append((180,split_lat_1))
            poly_2_pts = [(-180,split_lat_1)]
            poly_2_pts.extend(poly_bound_pts[split_ixs[0]+1:split_ixs[1]+1])
            poly_2_pts.append((-180,split_lat_2))
            poly_3_pts = [(180,split_lat_2)]
            poly_3_pts.extend(poly_bound_pts[split_ixs[1]+1:split_ixs[2]+1])
            poly_3_pts.append((180,split_lat_3))
            poly_2_pts.append((-180,split_lat_3))
            poly_2_pts.extend(poly_bound_pts[split_ixs[2]+1:split_ixs[3]+1])
            poly_2_pts.append((-180,split_lat_4))
            poly_1_pts.append((180,split_lat_4))
        
        # Case 3: EH bowl / WH triangles, descending orbit
        elif split_lat_4 < split_lat_1 < split_lat_2 < split_lat_3:
            poly_1_pts.append((-180,split_lat_1))
            poly_2_pts = [(180,split_lat_1)]
            poly_2_pts.extend(poly_bound_pts[split_ixs[0]+1:split_ixs[1]+1])
            poly_2_pts.append((180,split_lat_2))
            poly_3_pts = [(-180,split_lat_2)]
            poly_3_pts.extend(poly_bound_pts[split_ixs[1]+1:split_ixs[2]+1])
            poly_3_pts.append((-180,split_lat_3))
            poly_2_pts.append((180,split_lat_3))
            poly_2_pts.extend(poly_bound_pts[split_ixs[2]+1:split_ixs[3]+1])
            poly_2_pts.append((180,split_lat_4))
            poly_1_pts.append((-180,split_lat_4))
        
        # Case 4: WH bowl / EH triangles, descending orbit
        elif split_lat_4 < split_lat_3 < split_lat_2 < split_lat_1:
            poly_1_pts.append((-180,split_lat_1))
            poly_2_pts = [(180,split_lat_1)]
            poly_2_pts.extend(poly_bound_pts[split_ixs[0]+1:split_ixs[1]+1])
            poly_2_pts.append((180,split_lat_2))
            poly_1_pts.append((-180,split_lat_2))
            poly_1_pts.extend(poly_bound_pts[split_ixs[1]+1:split_ixs[2]+1])
            poly_1_pts.append((-180,split_lat_3))
            poly_3_pts = [(180,split_lat_3)]
            poly_3_pts.extend(poly_bound_pts[split_ixs[2]+1:split_ixs[3]+1])
            poly_3_pts.append((180,split_lat_4))
            poly_1_pts.append((-180,split_lat_4))
        
        poly_1_pts.extend(poly_bound_pts[split_ixs[3]+1:])
    
    split_polys_list = [poly_1_pts, poly_2_pts, poly_3_pts]
    
    return split_polys_list    


def TIRS_L1B_time_array_to_dt(TIRS_time_array):
    """
    Given a single TIRS observation time from the 'times_UTC' netCDF variable,
    represented as an array with shape (7,), return a python datetime object. 

    Parameters
    ----------
    TIRS_time_array : numpy.ndarray
        TIRS observation time represented as as an array with shape (7,). Array
        entries are: year, month, day, hour, minute, second, microsecond.

    Returns
    -------
    TIRS_dt : datetime.datetime
        TIRS observation time as python datetime.

    """
    TIRS_dt = datetime.datetime(
        int(TIRS_time_array[0]),
        int(TIRS_time_array[1]),
        int(TIRS_time_array[2]),
        int(TIRS_time_array[3]),
        int(TIRS_time_array[4]),
        int(TIRS_time_array[5]),
        int(TIRS_time_array[6])
        )
    
    return TIRS_dt