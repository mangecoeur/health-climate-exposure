import rtree.index

from joblib import Memory

from pathlib import Path

import numpy as np
import xarray as xr


from numba import jit
from shapely.geometry import Polygon
from shapely.geometry.base import BaseGeometry

from util import GridLookupResults

CACHE_BASE = Path('~/Data/python/health-climate-exposure/').expanduser()


memory = Memory(cachedir=CACHE_BASE)

def save_latlon_arrays(weather_file):
    latest_weather = xr.open_dataset(weather_file)
    lon = np.array(latest_weather.longitude)
    lat = np.array(latest_weather.latitude)
    np.save('longitudes.npy', lon)
    np.save('latitudes.npy', lat)


@jit
def generate_polygon_points(longitude, latitude):
    n_lon = len(longitude)
    n_lat = len(latitude)

    lat_m, lon_m = np.meshgrid(longitude, latitude)
    n_polys = (n_lat - 1) * (n_lon - 1)

    out = np.zeros((n_polys, 4, 2))

    temp_pol = np.zeros((4, 2))

    stacked_latlon = np.dstack([lat_m, lon_m])
    counter = 0  # FIXME no seriously, you are smart enough not to need this :/

    for lon_idx in range(n_lon - 1):
        for lat_idx in range(n_lat - 1):
            r = stacked_latlon[lat_idx:lat_idx + 2, lon_idx:lon_idx + 2, :]
            # r 2x2x2 - 2 lon cols, 2 lat rows, 2 height covering stacked lon, lat mgrids
            # 4 points A B C D
            # FIXME: should be possible to do this using rshape...
            temp_pol[0, :] = r[0, 0, :]
            temp_pol[1, :] = r[0, 1, :]
            temp_pol[2, :] = r[1, 1, :]
            temp_pol[3, :] = r[1, 0, :]

            # out[lon_idx + (lon_idx + 1) * lat_idx, :, :] = r
            out[counter, :, :] = temp_pol

            counter += 1

    return out


@memory.cache
def weather_file_grid(weather_dataset):
    """
    data is list of 4 points of quad

    D ---- C
    |      |
    |      |
    A ---- B
    :param data:
    :return:
    """
    # longitude = np.load(Path('data') / 'longitudes.npy')
    # latitude = np.load(Path('data') / 'latitudes.npy')
    # weather_dataset = xr.open_dataset(weather_file)
    longitude = np.array(weather_dataset.longitude)
    latitude = np.array(weather_dataset.latitude)

    # Since we want a WGS84 grid, fold values >180 into range -180 to 0
    # longitude -= 180
    # To avoid making a mess when we then turn the points into grid squares need to re-order
    # the longitude so that it goes from min to max
    longitude[longitude > 180] -= 360
    # longitude = np.concatenate([longitude[longitude < 0], longitude[longitude >= 0]])
    longitude = np.sort(longitude)
    latitude = np.sort(latitude)

    return generate_polygon_points(longitude, latitude)

def build_save_index(rects):
    print('creating file index')
    idx = rtree.index.Index('era_interim')

    for r_id, rect in enumerate(rects):
        point = rect[0, :]

        lo_shift = (rect[1, 0] - rect[0, 0]) / 2
        la_shift = (rect[3, 1] - rect[0, 1]) / 2

        rect[:, 0] -= lo_shift
        rect[:, 1] -= la_shift

        bbox = (rect[0, 0], rect[0, 1], rect[2, 0], rect[2, 1])

        idx.insert(r_id, bbox, {'point': point, 'bbox': rect})
    return idx

def index_points_and_bbox():
    idx = rtree.index.Index('era_interim')
    return idx


_INDEX = None

def get_index():
    global _INDEX
    if not Path('era_interim.dat').is_file():
        rects = weather_file_grid()
        _INDEX = build_save_index(rects)
    else:
        _INDEX = index_points_and_bbox()
    return _INDEX


def find_shape_in_index(shape: BaseGeometry, index):
    # (minx, miny, maxx, maxy)
    bounds = shape.bounds

    lon = []
    lat = []
    weights = []

    for result in index.intersection(bounds, objects="raw"):
        point = result['point']
        rect = result['bbox']
        poly = Polygon([(r[0], r[1]) for r in rect])

        inters = shape.intersection(poly)
        # Intersection with bounds might not intersect with detailed poly
        if inters:
            weights.append(inters.area)
            lon.append(point[0])
            lat.append(point[1])

    lon = np.array(lon)
    lat = np.array(lat)
    weights = np.array(weights)
    weights = weights / np.sum(weights)
    return GridLookupResults(lon, lat, weights)


def find_in_era(shape):
    return find_shape_in_index(shape, get_index())

# import cartopy
#
#
# SHAPEFILES_SRC = DATA_SRC / 'GIS' / 'world'
#
# WEATHER_FILE = DATA_SRC / 'weather' / 'ecmwf' / 'ecmwf_era_interim_2014-2017.nc'
#
#
# COUNTRY_POLY_SHP = SHAPEFILES_SRC / 'ne_10m_admin_0_countries' / 'ne_10m_admin_0_countries.shp'
#
# print('Open files')
# countries = list(cartopy.io.shapereader.Reader(str(COUNTRY_POLY_SHP)).geometries())
#
# print('Era search')
# res = find_in_era(countries[10])
# print(res)