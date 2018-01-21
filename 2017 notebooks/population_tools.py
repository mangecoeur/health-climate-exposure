"""
Interpolate the NASA gridded populations, given every 5 years
to get population grid for every year. Do this with a simple 
linear interpolation.
"""
from contextlib import AbstractContextManager
from enum import Enum
from pathlib import Path

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
    out = np.roll(out, -width // 2, axis=1)

    # Create dataset
    # ds = xr.Dataset({'population': (['latitude', 'longitude', 'time'], out)},
    #                 coords={'longitude': pop_x,
    #                         'latitude': pop_y,
    #                         'time': pd.date_range('2000-01-01', periods=out.shape[2], freq='AS')
    #                         })

    ds = xr.Dataset({'population': (['latitude', 'longitude', 'year'], out)},
                    coords={'longitude': pop_x,
                            'latitude': pop_y,
                            'year': np.arange(2000, 2000 + out.shape[2], 1, dtype=np.int32)
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
        first = population[::2, :]
        second = population[1::2, :]
        if second.shape[0] < first.shape[0]:
            # missing a row, need to 'wrap'- just double the values from 'first' as an appoximation
            second = np.vstack((second, first[-1, :]))
        population = first + second
        # Sum every other column
        if second.shape[1] < first.shape[1]:
            # missing a row, need to 'wrap'- just double the values from 'first' as an appoximation
            second = np.hstack((second, first[:, -1]))
        first = population[:, ::2]
        second = population[:, 1::2]
        # population = population[:, ::2] + population[:, 1::2]
        population = first + second

        if how == 'mean':
            population = population / 4
        # Output affine scaled by 2
        trns = Affine(trns.a * 2, trns.b, trns.c, trns.d, trns.e * 2, trns.f)
    # Reduction to 1/4 of the original size already makes life much easier
    save_population_geotiff(trns, pop_meta, population, year)


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

        # derez_population(population_file, year, 4, how)

            executor.submit(derez_population, population_file, year, 4, how)


@jit(cache=True)
def reproject_to(shape, src_data, src_affine, out_affine, crs, out_lat, out_lon):
    reproj = np.empty(shape=shape[0:2])
    # TODO do you really need rasterio? could just use scipy interp.
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


def get_water_mask(target, file_path):
    """
    Get the water mask on 0-360 lon range

    Args:
        target:

    Returns:

    """

    with rasterio.open(str(file_path)) as pop:
        pop_mask = pop.read(1)


        new_mask = np.empty(shape=(len(target.latitude),
                                   len(target.longitude)))

        new_aff = get_affine(target)
        # Override xform to ensure -180 to 180 range
        transform = Affine(new_aff.a, 0, -180, 0, new_aff.e, new_aff.f)
        reproject(
            pop_mask, new_mask,
            src_transform=pop.transform,
            dst_transform=transform,
            src_crs=pop.crs,
            dst_crs=pop.crs,
            resample=Resampling.bilinear)

    # Roll to fix -180 to 180 vs 0 to 360 convention

    width = new_mask.shape[1]
    new_mask = np.roll(new_mask, -width // 2, axis=1)
    # TODO should this be as int or should we set 0 to nan
    new_mask[new_mask == 0] = np.nan
    return new_mask


class PopulationType(Enum):
    count = POP_DATA_SRC / 'population_count_2000-2020.nc'
    density = POP_DATA_SRC / 'population_density_2000-2020.nc'


DEFAULT_FILE = POP_DATA_SRC / 'population_count_2000-2020.nc'
# DEFAULT_FILE = POP_DATA_SRC / 'population_count_2000-2020_highres.nc'



class PopulationProjector(AbstractContextManager):
    def __init__(self, population_file,
                 water_mask_file=None,
                 mask_empty=True):
        self.crs = CRS({'init': 'epsg:4326'})

        # TODO don't force add this config path....
        pop_file = POP_DATA_SRC / population_file
        # self.data: xr.Dataset = xr.open_dataset(str(pop_file), chunks={'time': 2})
        self.data: xr.Dataset = xr.open_dataarray(str(pop_file), chunks={'year': 2})

        if 'lon' in self.data.dims:
            self.data = self.data.rename({'lon': 'longitude', 'lat': 'latitude'})

        # self.data['time'] =  self.data['time.year']
        # self.data =  self.data.rename({'time': 'year'})
        self.affine = get_affine(self.data)

        # water_mask_path = POP_DATA_SRC / 'water_mask_eightres.tif'
        if water_mask_file is None:
            water_mask_path = POP_DATA_SRC / 'water_mask_sixteenres.tif'
        else:
            water_mask_path = POP_DATA_SRC / water_mask_file

        self.water_mask = get_water_mask(self.data, water_mask_path)
        # self.water_mask.shape = (*self.water_mask.shape, 1)
        self.mask_empty = mask_empty

        if self.mask_empty:
            self.data = self.data_empty_masked

    @property
    def data_water_masked(self) -> xr.Dataset:
        """
        Returns:
            Population with areas of water and ice replaced with NaN
        """

        return self.data * self.water_mask

    @property
    def data_empty_masked(self) -> xr.Dataset:
        """
        Returns:
            Population with areas of water and ice replaced with NaN
        """
        da = xr.DataArray(self.water_mask, coords=[self.data.latitude, self.data.longitude],
                          dims=['latitude', 'longitude'],
                          name='water_mask')
        return da * self.data.where(self.data > 1e-08)

    def rasterize_data(self, table, key, affine=None):
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
            table: GeoPandas table
            key: Table column name

        Returns:

        """
        # NOTE: it seems that rasterize wants -180 to 180 :/
        # Figure out the transform from the geopandas to the raster
        # affine = Affine(self.population_affine.a, 0, px, 0, self.population_affine.e, self.population_affine.f)
        affine = affine if affine else self.affine

        raster = features.rasterize(
            ((r.geometry, r[key]) for _, r in table.iterrows()),
            out_shape=self.data.shape[:2],
            transform=affine
        )
        # Roll the result to fix affine oddity
        raster = np.roll(raster, -raster.shape[1] // 2, axis=1)

        return raster

    def project(self, year, param: xr.DataArray):
        """
        Project param onto the population grid and multiply by population values
        
        Args:
            year: 
            param: 

        Returns:

        """

        lon = param.longitude if 'longitude' in param.coords else param.lon
        lat = param.latitude if 'latitude' in param.coords else param.lat
        dx = (lon[1] - lon[0]).values
        px = lon[0].values

        dy = (lat[1] - lat[0]).values
        py = lat[0].values

        input_affine = Affine(dx, 0, px, 0, dy, py)

        if self.mask_empty:
            # projected = projected * self.water_mask
            pop_year = self.data_empty_masked.sel(year=year)
        else:
            pop_year = self.data.sel(year=year)

        projected = reproject_to(pop_year.shape, param,
                                 input_affine, self.affine, self.crs,
                                 pop_year.latitude, pop_year.longitude)

        projected = pop_year * projected

        return projected

    def project_param(self, param: xr.DataArray):
        """
        Project parameter to the population grid, but don't multiply
        
        Args:
            year: 
            param: 

        Returns:

        """

        lon = param.longitude if 'longitude' in param.coords else param.lon
        lat = param.latitude if 'latitude' in param.coords else param.lat
        dx = (lon[1] - lon[0]).values
        px = lon[0].values

        dy = (lat[1] - lat[0]).values
        py = lat[0].values

        input_affine = Affine(dx, 0, px, 0, dy, py)
        return reproject_to(self.data.shape, param, input_affine, self.affine, self.crs,
                            self.data.latitude, self.data.longitude)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.data.close()


# REZ_FIX = 'eightres'
REZ_FIX = 'sixteenres'

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
