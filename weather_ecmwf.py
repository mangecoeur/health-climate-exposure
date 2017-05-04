from collections import namedtuple
from enum import Enum
from pathlib import Path

import numpy as np
import pandas as pd
import rtree
import xarray as xr
from sqlalchemy import text

import spatial_lookup
from util import GridLookupResults


def _weighted_grid_lookup(code, sql, db):
    geo_data = db.execute(text(sql), code=code).fetchall()
    if len(geo_data) == 0:
        return GridLookupResults([], [], [])
    lon, lat, weights = tuple(zip(*geo_data))
    lon = np.array(lon)
    lat = np.array(lat)
    weights = np.array(weights) / np.sum(weights)

    return GridLookupResults(lon, lat, weights)


def lookup_country_grid(iso3_country_code, db, grid_table='era_interim_grid') -> GridLookupResults:
    country_sql = """
    SELECT st_X(grid_points), st_y(grid_points), 
           (st_area(st_intersection(g.grid_polygons, p.geom)::GEOGRAPHY)/st_area(p.geom::GEOGRAPHY)) AS pct_area
    FROM gis.{grid_table} AS g, gis.natural_earth_countries AS p
    WHERE st_intersects(g.grid_polygons, p.geom)
    AND p.iso_a3 =:code;
    """.format(grid_table=grid_table)

    return _weighted_grid_lookup(iso3_country_code, country_sql, db)


_Row = namedtuple('_Row', ['ncdf_name', 'common_name'])

# Manually extracted from the weather files.
# python compat var long name, netcdf file name, common name
WEATHER_NAMES = [
    _Row('t2m', 'temperature_2m'),
    _Row('tp', 'precipitation'),
    _Row('d2m', 'temperature_dewpoint'),
    _Row('sp', 'surface_pressure'),
    _Row('2T_GDS4_SFC', 'temperature_2m'),
    _Row('g4_lat_1', 'latitude'),
    _Row('g4_lon_2', 'longitude'),
    _Row('initial_time0_hours', 'time')
]


def weather_dataset(root_path, rename=True):
    data = xr.open_dataset(str(root_path))
    if rename:
        data = data.rename({r.ncdf_name: r.common_name for r in WEATHER_NAMES if r.ncdf_name in data.variables})
    return data


def climatology_dataset(file_path, rename=True, decode_time=True) -> xr.Dataset:
    climatology = xr.open_dataset(str(file_path), decode_times=decode_time)

    if 'initial_time0_encoded' in climatology.variables:
        synthetic_time = pd.date_range(start='1999-01-01', periods=len(climatology.initial_time0_hours.values),
                                       freq='6H')
        climatology['initial_time0_hours'] = synthetic_time
        climatology = climatology.drop(['initial_time0_encoded', 'initial_time0'])

    if rename:
        climatology = climatology.rename({r.ncdf_name: r.common_name
                                          for r in WEATHER_NAMES if r.ncdf_name in climatology.variables})

    return climatology


def weighted_regional_timeseries(dataset: xr.Dataset, lon, lat, weights=None, start_timestamp=None,
                                 end_timestamp=None) -> pd.DataFrame:
    """

    :param dataset:
    :param start_timestamp:
    :param end_timestamp:
    :param lon: lon and lats should be the same size and follow WGS84 conventions
    :param lat: lat
    :param weights: optional, defaults to equal weighting
    :param variables: optionally override the default names as found in the netCDF file

    :return:
       Dataframe timeseries of weather data averaged over the series selected by lon/lat/time indexers
       if weights are given the result is a weighted average

    """

    if len(lon) != len(lat):
        raise ValueError('Both lon and lat must tbe the same size')


    # ECMWF netCDF file has longitude 0-360 instead of -180-180 as in SRID 4326
    # So convert negative lon to range 0-360
    lon = np.asarray(lon)
    if np.min(lon) < 0:
        lon[lon < 0] = 360 + lon[lon < 0]
    lat = np.asarray(lat)

    # Make sure weights has the right shape to be broadcast against the dataset
    if weights is None:
        weights = np.ones_like(lat) / len(lat)

    weights = np.array(weights)

    # Note: because the geographic area calculations aren't always super accurate we need to allow
    # some tolerance for the weights - only good to 4 decimal places
    if not np.isclose(np.sum(weights), 1.0, 0.0001):
        raise ValueError('the given weights do not add up to 1±0.0001')


    # -- THIS IS WHERE THE SELECTION MAGIC HAPPENS --
    if start_timestamp or end_timestamp:
        time_sel = slice(start_timestamp, end_timestamp)
        selection = dataset.sel(time=time_sel)
    else:
        selection = dataset

    selection = selection.sel_points(method='nearest',
                                     tolerance=0.1,
                                     longitude=lon,
                                     latitude=lat).set_coords(['longitude', 'latitude', 'time'])
    # --

    if len(selection.points) != len(weights):
        raise ValueError('The number of supplied weights dont match the number of selected grid points')

    # Turn weights into a labeled dataarray so that it automatically gets broadcast
    # along the 'points' dimension
    weights = xr.DataArray(weights, [('points', selection['points'].values)])
    selection *= weights
    return selection.sum(dim='points')


class GridTable(Enum):
    era_interim_grid = 'era_interim_grid'
    era_climatology_grid = 'era_climatology_grid'


def lookup_country_data(db, dataset, country_iso_3, start_date=None, end_date=None,
                        resample=None, grid_table=GridTable.era_interim_grid):
    """
    
    Args:
        db: 
        dataset: 
        country_iso_3: 
        start_date: 
        end_date: 
        resample: 
        grid_table: 

    Returns:

    """
    grid_table = GridTable(grid_table)
    geo_lookup = lookup_country_grid(country_iso_3, db, grid_table.name)

    data = weighted_regional_timeseries(dataset, lon=geo_lookup.lon, lat=geo_lookup.lat, weights=geo_lookup.weights,
                                        start_timestamp=start_date, end_timestamp=end_date)

    if resample is not None:
        data = data.resample(resample, dim='time')

    return data


def lookup_shape_data(dataset, shape, start_date=None, end_date=None, resample=None):
    """
    
    Args:
        dataset: 
        shape: 
        start_date: 
        end_date: 
        resample: 

    Returns:

    """
    index = get_index(dataset)
    geo_lookup = spatial_lookup.find_shape_in_index(shape, index)
    data = weighted_regional_timeseries(dataset, lon=geo_lookup.lon, lat=geo_lookup.lat, weights=geo_lookup.weights,
                                        start_timestamp=start_date, end_timestamp=end_date)

    if resample is not None:
        data = data.resample(resample, dim='time')

    return data


_INDEX = None

def get_index(dataset):
    # TODO per-dataset index
    global _INDEX
    if not Path('era_interim.dat').is_file():
        rects = spatial_lookup.weather_file_grid(dataset)
        _INDEX = spatial_lookup.build_save_index(rects)
    else:
        _INDEX = rtree.index.Index('era_interim')
    return _INDEX