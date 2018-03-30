from pathlib import Path
from datetime import date

import numpy as np
import scipy as sp
import pandas as pd
import xarray as xr

from numba import jit, float64, float32

DATA_SRC = Path('~/Data/').expanduser()
WEATHER_SRC = DATA_SRC / 'weather'
OUT_FOLDER = WEATHER_SRC / 'ecmwf' / 'daily_wbgt'

from collections import namedtuple

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


def open_mfdataset(path_pattern, rename=True, **kwargs):
    data = xr.open_mfdataset(str(path_pattern), engine='scipy', **kwargs)
    if rename:
        data = data.rename({r.ncdf_name: r.common_name for r in WEATHER_NAMES if r.ncdf_name in data.variables})
    return data





@jit(nopython=True, nogil=True)
def calculate_relative_humidity(temperature, temperature_dewpoint):
    """
    RH = 100% x (E/Es)

    where, according to an approximation of the Clausius-Clapeyron equation:
    
    E = E0 x exp[(L/Rv) x {(1/T0) - (1/Td)}] and
    
    Es = E0 x exp[(L/Rv) x {(1/T0) - (1/T)}]
    
    where E0 = 0.611 kPa, (L/Rv) = 5423 K (in Kelvin, over a flat surface of water), T0 = 273 K (Kelvin)

    Returns:

    """

    E0 = 0.611
    t_0 = 273.15
    # (L / Rv) = 5423
    L_over_RV = 5423
    E = E0 * np.exp((L_over_RV) * ((1 / t_0) - (1 / temperature_dewpoint)))

    Es = E0 * np.exp((L_over_RV) * ((1 / t_0) - (1 / temperature)))

    RH = 100 * (E / Es)
    return RH

@jit(nopython=True, nogil=True)
def calculate_wbt_t_dew(t_ref, temperature_dewpoint, surface_pressure):
    """Empirical calculation of wet bulb temperature from temperature, humidity, and pressure

    Args:
        t_ref: Dry bulb near-surface air temperature (K)
        relative_humidity: Relative humidity (%)
        surface_pressure: Surface air pressure (Pa)

    Returns:
        Wet Bulb Temperature
    """

    # Empirical formula for e_sat
    a = (-2991.2729 / t_ref ** 2)
    b = -6017.0128 / t_ref
    c = -0.028354721 * t_ref
    d = 1.7838301E-5 * t_ref ** 2
    e = -8.4150417E-10 * t_ref ** 3
    f = 4.4412543E-13 * t_ref ** 4
    g = 2.858487 * np.log(t_ref)
    e_sat = np.exp(18.87643854 + a + b + c + d + e + f + g) / 100
    surface_pressure = surface_pressure / 100
    w_sat = 621.97 * e_sat / (surface_pressure - e_sat)
    
    relative_humidity = calculate_relative_humidity(t_ref, temperature_dewpoint)
    
    humidity_frac = relative_humidity / 100
    w = humidity_frac * w_sat
    t_l = 1 / (1 / (t_ref - 55) - np.log(humidity_frac) / 2840) + 55
    t_e = t_ref * (1000 / surface_pressure) ** (0.2854 * (1 - 0.28E-3 * w)) * np.exp(
        (3.376 / t_l - 0.00254) * w * (1 + 0.81E-3 * w))

    wbt = 45.114 - 51.489 * (t_e / 273.15) ** (-3.504)  # in ˚C
    # Standardize on kelvin for sanity
    return wbt + 273.15


@jit(nopython=True, nogil=True)
def calculate_wbgt_t_dew(t_ref, temperature_dewpoint, surface_pressure):
    wbt = calculate_wbt_t_dew(t_ref, temperature_dewpoint, surface_pressure)
    wbgt = 0.7 * (wbt - 273.15) + 0.3 * (t_ref - 273.15)
    return wbgt + 273.15
    
    

def main():
    weather = open_mfdataset(WEATHER_SRC / 'ecmwf' / 'daily_test/19*.nc',
                                           chunks={'time':100}
                                          )


    wbgt = xr.apply_ufunc(calculate_wbgt_t_dew,
                         weather.temperature_2m,  weather.temperature_dewpoint, weather.surface_pressure,
                         output_dtypes=[np.dtype(weather.temperature_2m)],
                         dask='parallelized'
                        )

    wbgt.name = 'wbgt'
    wbgt = wbgt.to_dataset()

    years, datasets = zip(*wbgt.groupby('time.year'))
    paths = [OUT_FOLDER / f'wbgt_{y}.nc' for y in years]
    xr.save_mfdataset(datasets, paths)
    
if __name__ == '__main__':
    main()