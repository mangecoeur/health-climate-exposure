import numpy as np
import xarray as xr

from config import GIS_DB, CLIMATOLOGY_FILE
from util import DB, generate_polygon_points, postgis_geom


def save_latlon_arrays():
    weather = xr.open_dataset(CLIMATOLOGY_FILE, decode_times=False)
    lon = np.array(weather.g4_lon_2)
    lat = np.array(weather.g4_lat_1)
    np.save('data/climatology_longitudes.npy', lon)
    np.save('data/climatology_latitudes.npy', lat)


def build_insert_era_interim(quad):
    point, poly = postgis_geom(quad)
    sql = f"""
    INSERT INTO gis.era_climatology_grid (grid_points, grid_polygons)
    VALUES (ST_GeomFromText('{point}', 4326), ST_GeomFromText('{poly}', 4326))
    """
    return sql



def era_interim_grid():

    weather = xr.open_dataset(CLIMATOLOGY_FILE, decode_times=False)
    longitude = np.array(weather.g4_lon_2)
    latitude = np.array(weather.g4_lat_1)


    # Since we want a WGS84 grid, fold values >180 into range -180 to 0
    # longitude -= 180
    # To avoid making a mess when we then turn the points into grid squares need to re-order
    # the longitude so that it goes from min to max
    longitude[longitude > 180] -= 360
    # longitude = np.concatenate([longitude[longitude < 0], longitude[longitude >= 0]])
    longitude = np.sort(longitude)

    return generate_polygon_points(longitude, latitude)


def build_era_interim(db):
    poly_points = era_interim_grid()

    create_table_era_interim(db)

    # Parallel(n_jobs=4)(delayed(do_insert)(quad) for quad in poly_points)
    # with ProcessPoolExecutor() as executor:
    #     executor.map(do_insert, poly_points)

    # TODO: could run this in parallel - but had problems doing so, so don't bother for now.
    for quad in poly_points:
        sql = build_insert_era_interim(quad)
        db.execute(sql)

    db.execute(
        'CREATE INDEX IF NOT EXISTS sidx_era_interim_world_grid ON gis.era_climatology_grid USING GIST ( grid_polygons );')


def create_table_era_interim(db):
    db.execute("""CREATE TABLE IF NOT EXISTS
      gis.era_climatology_grid(
        id SERIAL PRIMARY KEY NOT NULL,
    grid_points GEOMETRY(POINT, 4326),
    grid_polygons GEOMETRY(POLYGON, 4326)
    );"""
               )

    db.execute("""TRUNCATE TABLE gis.era_climatology_grid
    RESTART IDENTITY
    RESTRICT;""")


if __name__ == '__main__':
    gis_db = DB(GIS_DB)
    build_era_interim(gis_db)


