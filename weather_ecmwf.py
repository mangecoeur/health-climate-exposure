from collections import namedtuple, Callable
from contextlib import AbstractContextManager
from enum import Enum

import numpy as np
import pandas as pd
import xarray as xr
from sqlalchemy import text

from config import WEATHER_FILE, WEATHER_START, WEATHER_END, CLIMATOLOGY_FILE

GridLookupResults = namedtuple('GridLookupResults', ['lon', 'lat', 'weights'])


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
    _Row('2T_GDS4_SFC', 'temperature_2m')
]


def weather_dataset(root_path):
    return xr.open_dataset(root_path)


def climatology_dataset(file_path):
    climatology = xr.open_dataset(file_path, decode_times=False)
    synthetic_time = pd.date_range(start='1999-01-01', periods=len(climatology.initial_time0_hours.values), freq='6H')
    climatology['initial_time0_hours'] = synthetic_time
    climatology = climatology.drop(['initial_time0_encoded', 'initial_time0'])
    return climatology.rename({'g4_lat_1': 'latitude', 'g4_lon_2': 'longitude', 'initial_time0_hours': 'time'})


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
                        resample=None, rename=True, grid_table=GridTable.era_interim_grid):
    """

    :param db:
    :param dataset:
    :param lookup_fn:
    :param country_iso_3:
    :param start_date:
    :param end_date:
    :param grid_table:
    :param variables:
    :return:
    """
    grid_table = GridTable(grid_table)
    geo_lookup = lookup_country_grid(country_iso_3, db, grid_table.name)

    data = weighted_regional_timeseries(dataset, lon=geo_lookup.lon, lat=geo_lookup.lat, weights=geo_lookup.weights,
                                        start_timestamp=start_date, end_timestamp=end_date)

    if rename:
        renamer = {r.ncdf_name: r.common_name for r in WEATHER_NAMES if r.ncdf_name in data.variables}
        data = data.rename(renamer)

    if resample is not None:
        # Using Xarray resample method rather than pandas dataframe method
        data = data.resample(resample, dim='time')

    return data


# TODO decide whether I bother keeping any of this
class EraWeather(AbstractContextManager, Callable):
    """
    Wraps up an Xarray dataset configured to access the weather data.
    Behaves like a Dataset through delegation and can be used in a context manager -
    this is useful since setting up a dataset can be slow, and needs to be cleaned up
    after loading to avoid bugs/leaks. Using it in a context manager means you
    open it once and automatically clean up afterwards.

    Is also callable to enable 'functional style' processing - passes through to
    `lookup_site_data` method
    """

    def __init__(self, db, dset_start=None, dset_end=None, resample=None):
        """

        :param db:
        :param lookup_fn:
        :param dset_start:
        :param dset_end:
        :param site_conf:
        :param variables: list of 'common names' of variables to load, so that you can reduce the number
        of variables needed if you want.
        """
        self.dset_start = dset_start if dset_start else WEATHER_START
        self.dset_end = dset_end if dset_end else WEATHER_END

        self.db = db

        self.resample = resample

        self.weather = weather_dataset(WEATHER_FILE)

    def __exit__(self, type, value, traceback):
        self.close()

    def close(self):
        self.weather.close()

    def __call__(self, country_iso, start_date=None, end_date=None):
        return self.lookup_site_data(country_iso, start_date, end_date)

    def weather_for_country(self, country_iso, start_date=None, end_date=None) -> pd.DataFrame:
        """

        :param site:
        :param start_date:
        :param end_date:
        :param resample: can override the setting for the WeatherData object
        TODO: not sure if that is a good idea
        :return:
        """
        if start_date is None or end_date is None:
            start_date, end_date = self.dset_start, self.dset_end

        if start_date is None or end_date is None:
            raise ValueError('Must have a supply start/end date or a site with a valid time_range')

        data = lookup_country_data(self.db, self.weather, country_iso, start_date, end_date, self.resample)

        return data.to_dataframe()

    def get_timeseries(self, start_dt, end_dt, lon, lat) -> xr.Dataset:
        return self.weather.sel(time=slice(start_dt, end_dt)).sel_points(method='nearest',
                                                                         tolerance=0.1,
                                                                         lon=lon,
                                                                         lat=lat)

    def get_single_timeseries(self, start_date, end_date, lon, lat) -> xr.Dataset:
        selection = self.weather.sel(time=slice(start_date, end_date)).sel_points(method='nearest',
                                                                                  tolerance=0.1,
                                                                                  lon=[lon], lat=[lat])

        selection = selection.rename({r.ncdf_name: r.common_name for r in WEATHER_NAMES})

        selection = selection.mean(dim='points')

        return selection

    def __getattr__(self, attr):
        """
        Delegation pattern: any attr or method calls will be automatically delegated to
        the weather :class:`xarray.Dataset` if they are not found on the CfsrWeather object.

        .. note:
            If you chose to use the hack to include dew-point temperature, it will note be loaded
            when using dataset methods like :meth:`xarray.Datset.sel`
        :param attr:
        :return:
        """
        return getattr(self.weather, attr)
