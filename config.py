from pathlib import Path

DATA_SRC = Path('~/Data/').expanduser()
HUMANS_SRC = DATA_SRC / 'humans'
POP_DATA_SRC = HUMANS_SRC / 'population'
SHAPEFILES_SRC = DATA_SRC / 'GIS' / 'world'

ERA_2015 = DATA_SRC / 'weather' / 'ecmwf' / 'ecmwf_era_interim_2014-2017_compress.nc'
ERA_MONTHLY_FILE = DATA_SRC / 'weather' / 'ecmwf' / 'era_interim_monthly_means.nc'

CLIMATOLOGY_FILE = DATA_SRC / 'weather' / 'ecmwf' / 'era_climatology.nc'
CLIMATOLOGY_FILE_RESAMP = DATA_SRC / 'weather' / 'ecmwf' / 'era_climatology_resamp.nc'
CLIMATOLOGY_FILE_MONTHLY = DATA_SRC / 'weather' / 'ecmwf' / 'era_climatology_monthly.nc'


COUNTRY_POLY_SHP = SHAPEFILES_SRC / 'ne_10m_admin_0_countries' / 'ne_10m_admin_0_countries.shp'
# population_grid_2015_file = (POP_DATA_SRC / 'nasa_grid' /
#                              'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-2015' /
#                              'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_2015.tif')
GIS_DB = {
    'driver': 'psycopg2',
    'host': 'localhost',
    'username': 'jonathanchambers',
    'database': 'buildings',
    'port': 5432,
    'schema': 'gis'
}