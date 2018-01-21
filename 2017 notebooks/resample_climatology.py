import concurrent.futures
import numpy as np

from numba import jit
from scipy import interpolate

import xarray as xr

import weather_ecmwf
from config import CLIMATOLOGY_FILE, DATA_SRC, CLIMATOLOGY_FILE_RESAMP, ERA_MONTHLY_FILE


def do_resample(all_met_data, grid_lat_eq_mesh, grid_lon_eq_mesh, plat, plon, out_shape):
    output_data = np.empty((out_shape))
    print('resampling each slice of data')
    # TODO: this could be parallel...
    for i in range(out_shape[2]):
        met_data = all_met_data[:, :, i]
        met_resamp = resample_slice(grid_lat_eq_mesh, grid_lon_eq_mesh, met_data, plat, plon)
        output_data[:, :, i] = met_resamp

    return output_data

@jit
def resample_slice(weather_slice, new_lon_mesh, new_lat_mesh, plon, plat):
    original_bspline = interpolate.RectBivariateSpline(plat, plon, weather_slice)
    resamp = original_bspline(new_lat_mesh.ravel(), new_lon_mesh.ravel(), grid=False)
    resamp = resamp.reshape(new_lon_mesh.shape)
    return resamp


@jit
def resample():
    era_climatology = weather_ecmwf.climatology_dataset(CLIMATOLOGY_FILE, decode_time=False)
    era_weather = weather_ecmwf.weather_dataset(ERA_MONTHLY_FILE)
    plon = era_climatology.longitude

    # Have to flip the latitude dim to make resample happy.
    plat = np.sort(era_climatology.latitude)

    lon_grid, lat_grid = np.meshgrid( era_weather.longitude.values, era_weather.latitude.values)

    out = np.empty((len(era_climatology.time),
                    len(era_weather.latitude), len(era_weather.longitude)),
                   dtype=float)

    out = xr.DataArray(out, coords={'time': era_climatology.time,
                                    'longitude': era_weather.longitude,
                                    'latitude': era_weather.latitude
                                    }, dims=['time', 'latitude', 'longitude'])
    print('Resampling')
    for i, t in enumerate(era_climatology.time):
        time_slice = era_climatology.temperature_2m.sel(time=t)
        resamp = resample_slice(time_slice, lon_grid, lat_grid, plon, plat)
        # Flip the result back to match the input
        resamp = np.flip(resamp, axis=0)
        out.loc[t, : ,:] = resamp


    ds = xr.Dataset({'temperature_2m': out},
                    coords={'longitude': out.longitude,
                            'latitude': out.latitude,
                            'time': out.time
                            }
                    )
    # Flip the latitude coords to match ERA interim
    # ds['temperature_2m'] = np.flip(ds['temperature_2m'], 1)

    print('Saving')
    ds.to_netcdf(str(DATA_SRC / 'weather' / 'ecmwf' / 'era_climatology_resamp.nc'),
                 encoding={'temperature_2m': {'zlib':True}}
                 )


def climatology_to_monthly():
    era_climatology = weather_ecmwf.climatology_dataset(CLIMATOLOGY_FILE_RESAMP)

    print('Resampling to month start')
    era_climatology_monthly = era_climatology.resample('MS', dim='time')
    print('Saving')
    era_climatology_monthly.to_netcdf(str(DATA_SRC / 'weather' / 'ecmwf' / 'era_climatology_monthly.nc'),
                                      encoding={'temperature_2m': {'zlib': True}})

if __name__ == '__main__':
    resample()
    climatology_to_monthly()
