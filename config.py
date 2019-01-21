from pathlib import Path

DATA_SRC = Path('~/Data/').expanduser()
WEATHER_SRC = DATA_SRC / 'weather'
HUMANS_SRC = DATA_SRC / 'lancet'
POP_DATA_SRC = HUMANS_SRC / 'population'
SHAPEFILES_SRC = DATA_SRC / 'GIS' / 'world'

# ERA_MONTHLY_FILE = DATA_SRC / 'weather' / 'ecmwf' / 'era_interim_monthly_means.nc'
# ERA_MONTHLY_FILE = DATA_SRC / 'weather' / 'ecmwf' / 'monthly_means'

CLIMATOLOGY_FILE = DATA_SRC / 'weather' / 'ecmwf' / 'era_climatology.nc'
CLIMATOLOGY_FILE_RESAMP = DATA_SRC / 'weather' / 'ecmwf' / 'era_climatology_resamp.nc'


COUNTRY_POLY_SHP = SHAPEFILES_SRC / 'ne_10m_admin_0_countries' / 'ne_10m_admin_0_countries.shp'

GIS_DB = {
    'driver': 'psycopg2',
    'host': 'localhost',
    'username': 'jonathanchambers',
    'database': 'buildings',
    'port': 5432,
    'schema': 'gis'
}