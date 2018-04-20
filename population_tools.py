"""
Interpolate the NASA gridded populations, given every 5 years
to get population grid for every year. Do this with a simple 
linear interpolation.
"""
from contextlib import AbstractContextManager
from pathlib import Path

import datetime
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
from tqdm import tnrange

from config import POP_DATA_SRC


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
    population_file_path_two = _POP_SRC / _POPULATION_PATH_TEMPLATE.format(year=year_two,
                                                                           resolution=REZ_FIX)
    # population_file_path_one = POP_DATA_SRC / 'nasa_grid' / 'count' / f'population_{year_one}_{REZ_FIX}.tif'
    # population_file_path_two = POP_DATA_SRC / 'nasa_grid' / 'count' / f'population_{year_two}_{REZ_FIX}.tif'

    with rasterio.open(str(population_file_path_one)) as pop:
        population_one = pop.read(1)

    with rasterio.open(str(population_file_path_two)) as pop:
        population_two = pop.read(1)

    gradient = (population_two - population_one) / interval
    return population_one, gradient


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


def interp_to_netcdf():
    interval = 5
    height, width = get_shape()
    common_crs, common_trns = get_crs_and_affine()
    dx, _, px, _, dy, py, _, _, _ = common_trns

    # New affine lon 0-360
    era_compat_affine = Affine(dx, 0, 0, 0, dy, py)
    out = create_timeseries(height, interval, width)

    print('Init netcdf')

    pop_x = np.arange(0, width)
    pop_y = np.arange(0, height)

    # change the coords so they match ERA lon from 0 to 360
    dx, _, px, _, dy, py, _, _, _ = era_compat_affine
    pop_x = pop_x * dx + px
    pop_y = pop_y * dy + py
    out = np.roll(out, -width // 2, axis=1)

    ds = xr.Dataset({'population': (['latitude', 'longitude', 'year'], out)},
                    coords={'longitude': pop_x,
                            'latitude': pop_y,
                            'year': np.arange(2000, 2000 + out.shape[2], 1, dtype=np.int32)
                            })

    print(ds)

    print('Saving')
    ds.to_netcdf(str(POP_DATA_SRC / 'population_{type}_2000-2020_{rezfix}.nc'.format(type=POP_TYPE, rezfix=REZ_FIX)),
                 format='NETCDF4',
                 encoding={'population': {
                     'zlib': True,
                 }}
                 )
    print('Done')


def derez_population(population_file_path, n_iters=1, how='sum'):
    with rasterio.open(str(population_file_path)) as pop:
        print(pop.meta)
        pop_meta = pop.meta
        trns = pop.transform
        population = pop.read(1, masked=True)

    population.fill_value = 0
    population = population.filled()

    for i in range(n_iters):
        # Sum every other row
        first = population[::2, :]
        second = population[1::2, :]
        if second.shape[0] < first.shape[0]:
            # missing a row, need to 'wrap'- just double the values from 'first' as an aproximation
            second = np.vstack((second, first[-1, :]))
        population = first + second
        # Sum every other column
        if second.shape[1] < first.shape[1]:
            # missing a row, need to 'wrap'- just double the values from 'first' as an aproximation
            second = np.hstack((second, first[:, -1]))
        first = population[:, ::2]
        second = population[:, 1::2]
        # population = population[:, ::2] + population[:, 1::2]
        population = first + second

        if how == 'mean':
            population = population / 4
        # Output affine scaled by 2
        trns = Affine(trns.a * 2, trns.b, trns.c, trns.d, trns.e * 2, trns.f)
    return population, pop_meta, trns


def save_population_geotiff(trns, pop_meta, population, year):
    print('Saving')
    print(population.shape)
    with rasterio.open(str(_POP_SRC / _POPULATION_PATH_TEMPLATE.format(year=year,
                                                                       resolution=REZ_FIX)),
                       'w',
                       driver='GTiff',
                       height=population.shape[0],
                       width=population.shape[1],
                       count=1,
                       dtype='float32',
                       crs=pop_meta['crs'],
                       transform=trns,
                       compress='lzw') as new_dataset:
        new_dataset.write(population, 1)


def derez_population_and_save_geotiff(population_file, year, n_iters, how):
    population, pop_meta, trns = derez_population(population_file, n_iters=n_iters, how=how)
    save_population_geotiff(trns, pop_meta, population, year)


def do_derez(how='sum'):
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=5) as executor:
        for year in [2000, 2005, 2010, 2015, 2020]:
            print(f'Proc for {year}')
            population_file = (_POP_SRC /
                               _POP_ORIGINAL_FOLDER_TMPL.format(type=POP_TYPE, year=year) /
                               (_POP_ORIGINAL_TMPL.format(type=POP_TYPE, year=year) + '.tif'))

            executor.submit(derez_population_and_save_geotiff, population_file, year, N_ITERS_DEREZ, how)


def get_affine_latlon(lat, lon):
    dx = (lon[1] - lon[0]).values
    px = lon[0].values
    dy = (lat[1] - lat[0]).values
    py = lat[0].values
    return Affine(dx, 0, px, 0, dy, py)


def get_affine(data):
    """
    Get an Affine transform from a data object with longitude and latitude attributes
    Args:
        data: 

    Returns:
        Affine
    """
    lon = data.longitude if 'longitude' in data.coords else data.lon
    lat = data.latitude if 'latitude' in data.coords else data.lat

    return get_affine_latlon(lat, lon)


def rasterize_data(target_dataset, table, key, affine=None) -> np.ndarray:
    """
    Rasterize a geopandas table with a geometry column
    onto the population grid, setting the shape regions to the
    value in the column name 'key'

    Returns a raster that should match the population bits,
    on grid 0 to 360.

    Note that need some fiddling for this, not certain how the
    projection is working but seems to want to assume -180 to 180 input/output
    so have to roll to get ERA compatible layout

    Args:
        target_dataset: dataset with lat, lon, and 2D grid data which will be used as the rasterization target
        table: GeoPandas table
        key: Table column name
        affine: Transformation from pixel coordinates of image to the
        coordinate system of the input shapes.

    Returns:
        Numpy array of rasterized table data rasterized onto same grid as
        population data

    """
    # NOTE: it seems that rasterize wants -180 to 180, even though it should accept arbitrary affine :/
    # Figure out the transform from the geopandas to the raster
    # affine = Affine(self.population_affine.a, 0, px, 0, self.population_affine.e, self.population_affine.f)
    affine = affine if affine else get_affine(target_dataset)

    raster = features.rasterize(
        ((r.geometry, r[key]) for _, r in table.iterrows()),
        out_shape=target_dataset.shape[:2],
        transform=affine
    )
    # Roll the result to fix affine oddity
    raster = np.roll(raster, -raster.shape[1] // 2, axis=1)

    return raster


def reproject_to(shape, src_data, src_affine, out_affine, crs, out_lat, out_lon):
    reproj = np.empty(shape=shape)

    reproject(
        src_data.values, reproj,
        src_transform=src_affine,
        dst_transform=out_affine,
        src_crs=crs,
        dst_crs=crs,
        resample=Resampling.cubic_spline)

    # Wrap as DataArray to make sure the coordinates line up
    reproj = xr.DataArray(reproj,
                          coords={'latitude': out_lat,
                                  'longitude': out_lon
                                  },
                          dims=['latitude', 'longitude'])
    return reproj


def project_param(target: xr.DataArray, param: xr.DataArray, crs=None):
    """
    Project parameter to the target grid, assuming all are in
    CRS 4326

    Args:
        target: Target DataArray with latitude and longitude coordinates
        param: 2D DataArray

    Returns:

    """
    target_affine = get_affine(target)
    input_affine = get_affine(param)

    crs = crs if crs else CRS({'init': 'epsg:4326'})

    # return reproject_to((1, *target.shape), param, input_affine, target_affine, crs,
    #                     target.latitude, target.longitude)
    return reproject_to(target.shape, param, input_affine, target_affine, crs,
                        target.latitude, target.longitude)


DEFAULT_POP_FILE = POP_DATA_SRC / 'population_count_2000-2020_eightres.nc'


def load_masked_population(population_file=DEFAULT_POP_FILE):
    """
    Shortcut to load the population file with the mask to remove areas of
    water and non-populated areas (population < 1e-8 to account for rouding errors)

    Args:
        population_file: path to population input file which includes 'population'
         and 'water_mask' variables

    Returns:
        population dataset including masked population data and water mask
    """
    data: xr.Dataset = xr.open_dataset(str(population_file),
                                       chunks={'year': 1})

    if 'lon' in data.dims:
        data = data.rename({'lon': 'longitude', 'lat': 'latitude'})

    data['empty_mask'] = data.water_mask * (data.population > 1e-08)
    data['population'] = data.water_mask * data.population.where(data.population > 1e-08)
    return data


def project_to_population(data: xr.DataArray,
                          weights=None, norm=False, start_year=2000, end_year=None,
                          population_file=DEFAULT_POP_FILE,
                          get_ts=True  # Shortcut, if you want the 3D grid rather then just the time series
                          ):
    """
    Resample a gridded timeseries dataset to match the grid of the given population file
    and multiply the values in each cell by the corresponding population, optionally
    normalising by the total population

    Args:
        data: the yearly data to project to the population grid
        weights: an optional grid with same dimensions as the population grid containing
        weightings for each cell, e.g. the fraction of the population to be considered for
        the given data
        norm (bool): whether to normalise by the total population, default False
        start_year:
        end_year:
        population_file: name of the population NetCDF file to use

    Returns:
        time series of population weighted and optionally `weights` weighted time series
    """
    # For convenience, accpet bot lat/lon and latitude/longitude
    if 'lon' in data.dims:
        data = data.rename({
            'lat': 'latitude',
            'lon': 'longitude'
        })


    data_time_dim = 'year'
    if data_time_dim not in data.dims:
        if 'time' not in data.dims:
            raise ValueError('data input must have a time-dimension named "time" or "year"')
        data_time_dim = 'time'

    if 'year' not in data.dims:
        data['year'] = data['time.year']

    with load_masked_population(population_file) as pop_data:

        if end_year is None:
            # Default to end of data
            end_year = np.asscalar(data['year'][-1])
            if 'time' in data.dims:
                end_time = data['time'][-1]


        if start_year is None:
            # If no start year is given, give the most recent of the data or the population data
            start_year = np.max(np.asscalar(data['year'][0]),
                                np.asscalar(pop_data['year'][0]))

        if weights is not None:
            pop_in = (pop_data.population * weights)
        else:
            pop_in = pop_data.population

        target_affine = get_affine(pop_data)
        input_affine = get_affine(data)

        def resample_ufunc(data_chunk, pop_chunk):
            out = np.empty(shape=pop_chunk.shape)
            crs = CRS({'init': 'epsg:4326'})

            for i in range(data_chunk.shape[0]):
                data_slice = data_chunk[i, :, :]
                out_slice = np.empty(shape=pop_chunk.shape[1:])
                reproject(
                    data_slice, out_slice,
                    src_transform=input_affine,
                    dst_transform=target_affine,
                    src_crs=crs,
                    dst_crs=crs,
                    resample=Resampling.cubic_spline)

                out[i, :, :] = out_slice * pop_chunk[i, :, :]

            return out

        # Select only the matching years so that apply_ufunc can join
        # on the year axis
        data_in = data.sel(year=slice(start_year, end_year))
        pop_in = pop_in.sel(year=slice(start_year, end_year))

        # Rename the lat/lon on the population so we can join the datasets on
        # year but make xarray understand that the lat/lon coords are different
        # for the different resolution datasets, otherwise it will try to join
        # using the coords but since we are resampling different resolutions
        # this will fail since they are different sizes
        pop_in = pop_in.rename({'latitude': 'pop_lat', 'longitude': 'pop_lon'})

        exposures = xr.apply_ufunc(resample_ufunc,
                                   data_in, pop_in,
                                   input_core_dims=[['latitude', 'longitude'], ['pop_lat', 'pop_lon']],
                                   output_core_dims=[['pop_lat', 'pop_lon']],
                                   dataset_join='outer',
                                   dask='parallelized',
                                   output_dtypes=[np.float64]
                                   )

        # We used the 'pop_lon/lat' as output dimensions to avoid name conflict, now rename to the normal name
        exposures = exposures.rename({'pop_lat': 'latitude',
                                      'pop_lon': 'longitude'})

        if norm:
            exposures = exposures / pop_in.sum(dim=['pop_lat', 'pop_lon'], skipna=True)

        if get_ts:
            exposures_ts = exposures.sum(dim=['latitude', 'longitude'],
                                         skipna=True)
            return exposures_ts

        else:
            return exposures


# def project_to_population(data, weights=None, norm=False, start_year=2000, end_year=None,
#                           population_file=DEFAULT_POP_FILE,
#                           get_ts=True  # Shortcut, if you want the 3D grid rather then just the time series
#                           ):
#     """
#     Resample a gridded timeseries dataset to match the grid of the given population file
#     and multiply the values in each cell by the corresponding population, optionally
#     normalising by the total population
#
#     Args:
#         data: the yearly data to project to the population grid
#         weights: an optional grid with same dimensions as the population grid containing
#         weightings for each cell, e.g. the fraction of the population to be considered for
#         the given data
#         norm (bool): whether to normalise by the total population, default False
#         start_year:
#         end_year:
#         population_file: name of the population NetCDF file to use
#
#     Returns:
#         time series of population weighted and optionally `weights` weighted time series
#     """
#     if end_year is None:
#         # Default to current year
#         end_year = datetime.datetime.now().year
#
#     with load_masked_population(population_file) as pop_data:
#         if weights is not None:
#             pop_sel = (pop_data.population * weights)
#         else:
#             pop_sel = pop_data.population
#         pop_sum = pop_sel.sum(dim=['latitude', 'longitude'], skipna=True)
#
#         target_shape = (len(pop_data.latitude), len(pop_data.longitude))
#         target_affine = get_affine(pop_data)
#         input_affine = get_affine(data)
#         crs = CRS({'init': 'epsg:4326'})
#
#         # This avoids trying to compute everything at once, which would use too much RAM
#         # TODO could re-write as ufunc to apply to pop_data, data pair and defer more stuff to play nice with Dask
#         # Idea is to make exposure constructing op return dask object
#         def do(year):
#             pop = pop_sel.sel(year=year).load()
#
#             proj = reproject_to(target_shape,
#                                 data.sel(year=year).load(),
#                                 input_affine, target_affine, crs,
#                                 pop_data.latitude, pop_data.longitude)
#
#             proj = proj * pop
#
#             return proj.compute().squeeze()
#
#         def do_ufunc(year):
#             pop = pop_sel.sel(year=year).load()
#
#             proj = reproject_to(target_shape,
#                                 data.sel(year=year).load(),
#                                 input_affine, target_affine, crs,
#                                 pop_data.latitude, pop_data.longitude)
#
#             proj = proj * pop
#
#             return proj.compute().squeeze()
#
#         xr.apply_ufunc(do_ufunc
#
#                        )
#         # TODO if we are only calculating the sum, do we need to bother creating a dataarray first.
#         exposures: xr.DataArray = xr.concat((do(year) for year in tnrange(start_year, end_year + 1)), dim='year')
#
#         if norm:
#             exposures = exposures / pop_sum
#
#         if get_ts:
#             exposures_ts = exposures.sum(dim=['latitude', 'longitude'],
#                                          skipna=True).compute()
#             return exposures_ts
#
#         else:
#             return exposures

# TODO clean this part up/move to notebook with explanations. Clean up global settings vars

N_ITERS_DEREZ = 3  # equivalent to 1/8th original resolution
REZ_FIX = 'eightres'
# REZ_FIX = 'sixteenres'

# POP_TYPE = 'density'
POP_TYPE = 'count'
_POPULATION_PATH_TEMPLATE = 'population_{year}_{resolution}.tif'
_POP_SRC = POP_DATA_SRC / 'nasa_grid' / POP_TYPE
# _POP_ORIGINAL_FOLDER_TMPL = 'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-{year}'
# _POP_ORIGINAL_TMPL = 'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-{year}'
_POP_ORIGINAL_FOLDER_TMPL = 'gpw-v4-population-{type}-adjusted-to-2015-unwpp-country-totals-{year}'
_POP_ORIGINAL_TMPL = 'gpw-v4-population-{type}-adjusted-to-2015-unwpp-country-totals_{year}'

if __name__ == '__main__':
    # do_derez(how='sum')
    interp_to_netcdf()
