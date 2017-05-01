"""
Interpolate the NASA gridded populations, given every 5 years
to get population grid for every year. Do this with a simple 
linear interpolation.
"""
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from affine import Affine
from numba import jit
from rasterio.enums import Resampling
from rasterio.warp import reproject

from config import POP_DATA_SRC

# REZ_FIX = 'quartres'
REZ_FIX = 'eightres'


def get_shape():
    year_one = 2000
    population_file_path_one = POP_DATA_SRC / 'nasa_grid' / 'count' / f'population_{year_one}_{REZ_FIX}.tif'

    with rasterio.open(str(population_file_path_one)) as pop:
        width = pop.width
        height = pop.height

    # rows, columns
    return height, width


def get_crs_and_affine():
    year_one = 2000
    population_file_path_one = POP_DATA_SRC / 'nasa_grid' / 'count' / f'population_{year_one}_{REZ_FIX}.tif'

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

    population_file_path_one = POP_DATA_SRC / 'nasa_grid' / 'count' / f'population_{year_one}_{REZ_FIX}.tif'
    population_file_path_two = POP_DATA_SRC / 'nasa_grid' / 'count' / f'population_{year_two}_{REZ_FIX}.tif'

    with rasterio.open(str(population_file_path_one)) as pop:
        population_one = pop.read(1)

    with rasterio.open(str(population_file_path_two)) as pop:
        population_two = pop.read(1)

    gradient = (population_two - population_one) / interval
    return population_one, gradient


def interp_to_gtiff(intercept, gradient, interval, year):
    """
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
    ds.to_netcdf(str(POP_DATA_SRC / 'population_2000-2020.nc'), format='NETCDF4',
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




def derez_population(population_file_path, year, n_iters=1):
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
        # Output affine scaled by 2
        trns = Affine(trns.a * 2, trns.b, trns.c, trns.d, trns.e * 2, trns.f)
    # Reduction to 1/4 of the original size already makes life much easier
    save_population_geotiff(trns, pop_meta, population, year)


def save_population_geotiff(newtrans, pop_meta, population, year):
    print('Saving')
    print(population.shape)
    with rasterio.open(str(POP_DATA_SRC / 'nasa_grid' / 'count' / f'population_{year}_{REZ_FIX}.tif'),
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


def do_derez():
    global year
    from concurrent.futures import ProcessPoolExecutor
    # year = 2005
    # pop_year = sys.argv[1]
    with ProcessPoolExecutor(max_workers=8) as executor:
        for year in [2000, 2005, 2010, 2015, 2020]:
            # for year in [2020]:
            print(f'Proc for {year}')

            population_file = (POP_DATA_SRC / 'nasa_grid' / 'count' /
                               f'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-{year}' /
                               f'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_{year}.tif')

            executor.submit(derez_population, population_file, year, 3)




            # derez_population(population_file, pop_year, 2)


if __name__ == '__main__':
    # do_derez()
    interp_to_netcdf()
