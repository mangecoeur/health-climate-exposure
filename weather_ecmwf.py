from collections import namedtuple

import pandas as pd
import xarray as xr

from ecmwfapi import ECMWFDataServer

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


def weather_mfdataset(root_path, rename=True):
    data = xr.open_mfdataset(str(root_path) + '/*.nc', engine='scipy')
    if rename:
        data = data.rename({r.ncdf_name: r.common_name for r in WEATHER_NAMES if r.ncdf_name in data.variables})
    return data




def any_weather_dataset(root_path, rename=True):
    data = xr.open_dataset(root_path)
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
