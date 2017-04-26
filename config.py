from pathlib import Path
from datetime import date

DATA_SRC = Path('~/Data/').expanduser()
POP_DATA_SRC = DATA_SRC / 'humans' / 'population'
SHAPEFILES_SRC = DATA_SRC / 'GIS' / 'world'

WEATHER_FILE = DATA_SRC / 'weather' / 'ecmwf' / 'ecmwf_era_interim_2014-2017.nc'
WEATHER_START = date(2015,1,1)
WEATHER_END = date(2017,1,1)

CLIMATOLOGY_FILE = DATA_SRC / 'weather' / 'ecmwf' / 'era_climatology.nc'


COUNTRY_POLY_SHP = SHAPEFILES_SRC / 'ne_10m_admin_0_countries' / 'ne_10m_admin_0_countries.shp'
population_grid_2015_file = (POP_DATA_SRC / 'nasa_grid' /
                             'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals-2015' /
                             'gpw-v4-population-count-adjusted-to-2015-unwpp-country-totals_2015.tif')
GIS_DB = {
    'driver': 'psycopg2',
    'host': 'localhost',
    'username': 'jonathanchambers',
    'database': 'buildings',
    'port': 5432,
    'schema': 'gis'
}