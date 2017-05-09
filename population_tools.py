"""
Interpolate the NASA gridded populations, given every 5 years
to get population grid for every year. Do this with a simple 
linear interpolation.
"""
from contextlib import AbstractContextManager
from enum import Enum

import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from affine import Affine
from numba import jit
from rasterio import features
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.warp import reproject

from config import POP_DATA_SRC

# REZ_FIX = 'quartres'
REZ_FIX = 'eightres'

POP_TYPE = 'density'
_POPULATION_PATH_TEMPLATE = 'population_{year}_{resolution}.tif'
# _POP_SRC = POP_DATA_SRC / 'nasa_grid' / 'count'
_POP_SRC = POP_DATA_SRC / 'nasa_grid' / POP_TYPE
# _POP_ORIGINAL_FOLDER_TMPL = 'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-{year}'
# _POP_ORIGINAL_TMPL = 'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-{year}'
_POP_ORIGINAL_FOLDER_TMPL = 'gpw-v4-population-{type}-adjusted-to-2015-unwpp-country-totals-{year}'
_POP_ORIGINAL_TMPL = 'gpw-v4-population-{type}-adjusted-to-2015-unwpp-country-totals_{year}'

def get_shape():
    year_one = 2000
    population_file_path_one = _POP_SRC / f'population_{year_one}_{REZ_FIX}.tif'

    with rasterio.open(str(population_file_path_one)) as pop:
        width = pop.width
        height = pop.height

    # rows, columns
    return height, width


def get_crs_and_affine():
    year_one = 2000
    population_file_path_one = _POP_SRC / f'population_{year_one}_{REZ_FIX}.tif'

    with rasterio.open(str(population_file_path_one)) as pop:
        pop_meta = pop.meta
        pop_trns = pop.transform

    return pop_meta['crs'], pop_trns


def get_era_compat_crs_affine():
    common_crs, common_trns = get_crs_and_affine()
    dx, _, px, _, dy, py, _, _, _ = common_trns

    # New affine lon 0-360
    era_compat_affine = Affine(dx, 0, 0, 0, dy, py)
    return common_crs, era_compat_affine

def lin_interp(year_one, interval):
    # Grids are at 5 year intervals
    # interval = 5
    year_two = year_one + interval

    population_file_path_one = _POP_SRC / _POPULATION_PATH_TEMPLATE.format(year=year_one,
                                                                           resolution=REZ_FIX)
    population_file_path_two =_POP_SRC / _POPULATION_PATH_TEMPLATE.format(year=year_two,
                                                                          resolution=REZ_FIX)
    # population_file_path_one = POP_DATA_SRC / 'nasa_grid' / 'count' / f'population_{year_one}_{REZ_FIX}.tif'
    # population_file_path_two = POP_DATA_SRC / 'nasa_grid' / 'count' / f'population_{year_two}_{REZ_FIX}.tif'

    with rasterio.open(str(population_file_path_one)) as pop:
        population_one = pop.read(1)

    with rasterio.open(str(population_file_path_two)) as pop:
        population_two = pop.read(1)

    gradient = (population_two - population_one) / interval
    return population_one, gradient



def interp_to_netcdf():
    interval = 5
    height, width = get_shape()
    common_crs, common_trns = get_crs_and_affine()
    dx, _, px, _, dy, py, _, _, _ = common_trns

    # New affine lon 0-360
    era_compat_affine = Affine(dx, 0, 0, 0, dy, py)
    out = create_timeseries(height, interval, width)

    print('Init netcdf')
    # Init the netcdf file
    pop_x = np.arange(0, width)
    pop_y = np.arange(0, height)

    # dx, _, px, _, dy, py, _, _, _ = common_trns

    # change the coords so they match ERA lon from 0 to 360
    dx, _, px, _, dy, py, _, _, _ = era_compat_affine
    pop_x = pop_x * dx + px
    pop_y = pop_y * dy + py
    out = np.roll(out, -width//2, axis=1)

    # Create dataset
    ds = xr.Dataset({'population': (['latitude', 'longitude', 'time'], out)},
                    coords={'longitude': pop_x,
                            'latitude': pop_y,
                            'time': pd.date_range('2000-01-01', periods=out.shape[2], freq='AS')
                            })

    print(ds)

    print('Saving')
    ds.to_netcdf(str(POP_DATA_SRC / 'population_{type}_2000-2020.nc'.format(type=POP_TYPE)), format='NETCDF4',
                 encoding={'population': {
                     'zlib': True,
                 }}
                 )
    print('Done')

@jit
def create_timeseries(height, interval, width):
    years = [2000, 2005, 2010, 2015]
    n_timesteps = len(years) * interval
    out = np.empty(shape=(height, width, n_timesteps))
    for yi, year in enumerate(years):
        print('Working on ', year)

        intercept, gradient = lin_interp(year, interval)

        for i in range(interval):
            out[:, :, yi * interval + i] = intercept + i * gradient
    return out




def derez_population(population_file_path, year, n_iters=1, how='sum'):
    with rasterio.open(str(population_file_path)) as pop:
        print(pop.meta)
        pop_meta = pop.meta
        trns = pop.transform
        population = pop.read(1, masked=True)
    population.fill_value = 0
    population = population.filled()

    for i in range(n_iters):
        # Sum every other row
        population = population[::2, :] + population[1::2, :]
        # Sum every other column
        population = population[:, ::2] + population[:, 1::2]
        if how == 'mean':
            population = population / 4
        # Output affine scaled by 2
        trns = Affine(trns.a * 2, trns.b, trns.c, trns.d, trns.e * 2, trns.f)
    # Reduction to 1/4 of the original size already makes life much easier
    save_population_geotiff(trns, pop_meta, population, year)


def save_population_geotiff(newtrans, pop_meta, population, year):
    print('Saving')
    print(population.shape)
    with rasterio.open(str(_POP_SRC / _POPULATION_PATH_TEMPLATE.format(year=year,
                                                                       resolution=REZ_FIX)),
                       'w',
                       driver='GTiff',
                       height=population.shape[0],
                       width=population.shape[1],
                       count=1,
                       dtype=population.dtype,
                       crs=pop_meta['crs'],
                       transform=newtrans,
                       compress='lzw') as new_dataset:
        new_dataset.write(population, 1)
        # new_dataset.write_mask(np.invert(population.mask))


def do_derez(how='sum'):
    global year
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=8) as executor:
        for year in [2000, 2005, 2010, 2015, 2020]:
            # for year in [2020]:
            print(f'Proc for {year}')
            population_file = (_POP_SRC /
                               _POP_ORIGINAL_FOLDER_TMPL.format(type=POP_TYPE, year=year) /
                               (_POP_ORIGINAL_TMPL.format(type=POP_TYPE, year=year) + '.tif'))

            executor.submit(derez_population, population_file, year, 3, how)

@jit
def reproject_to(shape, src_data, src_affine, out_affine, crs, out_lat, out_lon, out_time):
    reproj = np.empty(shape=shape[0:2])
    reproject(
        src_data.values, reproj,
        src_transform=src_affine,
        dst_transform=out_affine,
        src_crs=crs,
        dst_crs=crs,
        resample=Resampling.cubic_spline)

    # reproj.shape = (shape[0], shape[1], 1)
    # Wrap as DataArray to make sure the coordinates line up
    reproj = xr.DataArray(reproj,
                             coords={'latitude': out_lat,
                                     'longitude': out_lon,
                                     # 'time': out_time
                                     },
                             dims=['latitude', 'longitude',
                                   # 'time'
                                   ])
    return reproj


@jit
def project_to_population(year, src_data, pop_data, src_affine, out_affine, crs):
    pop_year = pop_data.population.sel(time='{year}'.format(year=year))

    resamp = reproject_to(pop_year.shape, src_data,
                          src_affine, out_affine, crs,
                          pop_year.latitude, pop_year.longitude, pop_year.time)
    res = pop_year * resamp
    return res


def get_affine(data):
    """
    Get an Affine transform from a data object with longitude and latitude attributes
    Args:
        data: 

    Returns:

    """
    lon = data.longitude
    lat = data.latitude
    dx = (lon[1] - lon[0]).values
    px = lon[0].values
    dy = (lat[1] - lat[0]).values
    py = lat[0].values
    return Affine(dx, 0, px, 0, dy, py)

#
# class PopulationType(Enum):
#     count = 'population_count_2000-2020.nc'
#     density = 'population_density_2000-2020.nc'
#

class PopulationProjector(AbstractContextManager):
    def __init__(self, population_file='population_count_2000-2020.nc', overlay_data=None, overlay_key='values'):
        self.crs = CRS({'init': 'epsg:4326'})

        pop_file = POP_DATA_SRC / population_file
        self.data = xr.open_dataset(str(pop_file))

        self.population_affine = get_affine(self.data)

        self.overlay = None

        # TODO should i rasterize an input here, or just pre-process it for each indicator
        if overlay_data:
            self.overlay = self.rasterize_data(overlay_data, overlay_key)


    def rasterize_data(self, table, key):
        raster = features.rasterize(
            ((r.geometry, r[key]) for _, r in table.iterrows()),
            out_shape=self.data.population.shape[:2],
            transform=self.population_affine
        )
        return raster

    def project(self, year, param: xr.DataArray):
        # TODO: might need to modify for when you also want to project demographic data

        lon = param.longitude if 'longitude' in param.coords else param.lon
        lat = param.latitude if 'latitude' in param.coords else param.lat
        dx = (lon[1] -lon[0]).values
        px = lon[0].values

        dy = (lat[1] - lat[0]).values
        py = lat[0].values

        input_affine = Affine(dx, 0, px, 0, dy, py)
        projected = project_to_population(year, param, self.data, input_affine, self.population_affine, self.crs)
        if self.overlay:
            projected = projected * self.overlay

        return projected

    def project_intermediate(self, year, param: xr.DataArray):
        # TODO: might need to modify for when you also want to project demographic data

        lon = param.longitude if 'longitude' in param.coords else param.lon
        lat = param.latitude if 'latitude' in param.coords else param.lat
        dx = (lon[1] -lon[0]).values
        px = lon[0].values

        dy = (lat[1] - lat[0]).values
        py = lat[0].values

        input_affine = Affine(dx, 0, px, 0, dy, py)
        return reproject_to(self.data.population.shape, param, input_affine, self.population_affine, self.crs,
                            self.data.latitude, self.data.longitude, [year])

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.data.close()


if __name__ == '__main__':
    do_derez(how='mean')
    interp_to_netcdf()


# --------------------
# Deprioritised work
# --------------------

def interp_to_gtiff(intercept, gradient, interval, year):
    """
    NOTE: preferred to do this to netcdf.
    
    Do the interpolation directly into the file to save on memory.
    
    Args:
        intercept: 
        gradient: 
        interval: 

    Returns:

    """
    common_crs, common_trns = get_crs_and_affine()

    with rasterio.open(str(POP_DATA_SRC / f'population_interp_{year}_{year+interval}.tif'),
                       'w',
                       driver='GTiff',
                       height=intercept.shape[0],
                       width=intercept.shape[1],
                       count=interval,
                       dtype=intercept.dtype,
                       crs=common_crs,
                       transform=common_trns,
                       compress='lzw') as new_dataset:
        for i in range(interval):
            pop_for_year = intercept + i * gradient
            new_dataset.write(pop_for_year, i + 1)  # gtiff bands numbered from 1

