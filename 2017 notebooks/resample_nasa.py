# coding: utf-8

# # Project NASA T anomaly onto population

# In[1]:

import numpy as np
import scipy as sp
import pandas as pd
import xarray as xr
import rasterio
from rasterio.warp import reproject, Resampling
from affine import Affine

from config import DATA_SRC, POP_DATA_SRC

nasa_giss_anom = DATA_SRC / 'weather' / 'nasa_giss' / 'air.2x2.1200.mon.anom.comb.nc'

nasa_trns = Affine(2, 0.0, -179, 0.0, -2, 89)


# TODO for intermediate years want to interp pop.
# ... or just use by country.
# TODO tweak to also support using summer months only
def dt_whole_year(year):
    population_file_path = POP_DATA_SRC / 'nasa_grid' / 'count' / f'population_{year}_quartres.tif'

    nasa_giss = xr.open_dataset(str(nasa_giss_anom))

    t_anomaly = nasa_giss.sel(time=f'{year}').air.mean(dim='time').values.astype('float32')

    with rasterio.open(str(population_file_path)) as pop:
        print(pop.meta)
        pop_meta = pop.meta
        pop_trns = pop.transform
        population = pop.read(1)

    common_crs = pop_meta['crs']


    newarr = np.empty(shape=population.shape)
    print('Reprojecting')

    reproject(
        t_anomaly, newarr,
        src_transform=nasa_trns,
        dst_transform=pop_trns,
        src_crs=common_crs,
        dst_crs=common_crs,
        resample=Resampling.bilinear)

    indicator = population * newarr

    with rasterio.open(str(DATA_SRC / 'lancet' / f'nasa_dt_indicator_{year}.tif'),
                       'w',
                       driver='GTiff',
                       height=indicator.shape[0],
                       width=indicator.shape[1],
                       count=1,
                       dtype=indicator.dtype,
                       crs=common_crs,
                       transform=pop_trns,
                       compress='lzw') as new_dataset:
        new_dataset.write(indicator, 1)


# year for population
# year = 2015
#
from concurrent.futures import ProcessPoolExecutor

with ProcessPoolExecutor(max_workers=8) as executor:
    executor.map(dt_whole_year, [2005,2010,2015])

# for year in [2005, 2010, 2015]:
    # dt_whole_year(year)
