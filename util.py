import os
import urllib.parse

import numpy as np
import pandas
import sqlalchemy
from numba import jit
from sqlalchemy import create_engine, MetaData


def _db_url(conf):
    # scheme://netloc/path;parameters?query#fragment
    u = conf['username'] if conf.get('username') else ''
    pwd = ':' + conf['password'] if conf.get('password') else ''

    p = ':' + str(conf['port']) if conf.get('port') else ''
    netloc = u + pwd + '@' + conf['host'] + p
    return urllib.parse.urlunparse(('postgres'+ '+' + conf.get('driver', 'psycopg2'),
                                    netloc,
                                    conf['database'],
                                    '', '', ''))


def create_engine_safe(db_url, schema):
    """
    Create a modified version of SqlAlchemy Engine which automatically
    scopes the connection to the given schema (so you don't need to prefix tables for that schema)
    and which automatically invalidates connections if the process ID (pid) changes
    such as when you use the multiprocessing module to perform DB access in
    parallel

    :param db_url:
    :param schema:
    :return:
    """
    # TODO: check whether schema is supported by DB first!
    engine = create_engine(db_url)
    meta = MetaData(db_url, schema=schema)
    meta.bind = engine

    @sqlalchemy.event.listens_for(engine, "connect")
    def connect(dbapi_connection, connection_record):
        connection_record.info['pid'] = os.getpid()
        # dbapi_connection.cursor().execute('set search_path={schema_name}, public'.format(schema_name=schema))

    @sqlalchemy.event.listens_for(engine, "checkout")
    def checkout(dbapi_connection, connection_record, connection_proxy):
        # Force the connection to use the schema search path.
        dbapi_connection.cursor().execute('set search_path={schema_name}, public'.format(schema_name=schema))

        pid = os.getpid()
        if connection_record.info['pid'] != pid:
            connection_record.connection = connection_proxy.connection = None
            raise sqlalchemy.exc.DisconnectionError(
                "Connection record belongs to pid %s, "
                "attempting to check out in pid %s" %
                (connection_record.info['pid'], pid)
            )

    # @event.listens_for(Pool, 'connect')
    # def set_search_path(dbapi_connection, conn_proxy):
    #     dbapi_connection.cursor().execute('set search_path={schema_name}, public'.format(schema_name=schema))

    return engine, meta


class DB(pandas.io.sql.SQLDatabase):
    def __init__(self, conf, **kwargs):
        """
        Class to help reading energy data from the Postgres database.
        Extends the PgSQLDatabase class, which itself is a subclass of the Pandas
        SQL handler. This setup uses postgres-specific performance optimisations
        (such as COPY FROM to load data) as well as some behaviour tweaks.
        Other utility functions exist to get, alter, and drop tables, add primary
        keys to existing tables, amonth others.

        This class initialises its own engine rather than requiring
        an engine object. It does this based on configuration information from
        the supplied config object.

        The database created is configured (in theory) such that connection
        pools are invalidated when a process is forked, so that if you use it
        in a multi-processing environment there shouldn't be any issues with
        trying to access pooled connections accros process boundaries.



        :param config: Config object
        """
        self.conf = conf
        schema = conf['schema'] if conf.get('schema') else kwargs.get('schema', 'public')

        db_url = _db_url(conf)
        self.engine, self.meta = create_engine_safe(db_url, schema)

        super().__init__(self.engine, meta=self.meta, **kwargs)
        self.schema = schema


@jit
def generate_polygon_points(longitude, latitude):
    n_lon = len(longitude)
    n_lat = len(latitude)

    lat_m, lon_m = np.meshgrid(longitude, latitude)
    n_polys = (n_lat - 1) * (n_lon - 1)

    out = np.zeros((n_polys, 4, 2))

    temp_pol = np.zeros((4, 2))

    stacked_latlon = np.dstack([lat_m, lon_m])

    counter = 0
    for lon_idx in range(n_lon - 1):
        for lat_idx in range(n_lat - 1):
            r = stacked_latlon[lat_idx:lat_idx + 2, lon_idx:lon_idx + 2, :]
            # r 2x2x2 - 2 lon cols, 2 lat rows, 2 height covering stacked lon, lat mgrids
            # 4 points A B C D
            # FIXME: should be possible to do this using rshape...
            temp_pol[0, :] = r[0, 0, :]
            temp_pol[1, :] = r[0, 1, :]
            temp_pol[2, :] = r[1, 1, :]
            temp_pol[3, :] = r[1, 0, :]

            # out[lon_idx + (lon_idx + 1) * lat_idx, :, :] = r
            out[counter, :, :] = temp_pol

            counter += 1

    return out


def postgis_rect(data):
    """
    data is list of 4 points of quad

    D ---- C
    |      |
    |      |
    A ---- B
    :param data:
    :return:
    """
    # TODO: shift all the points by half the width, half the height
    # Shift lon by (B-A)/2, shift lat by (D-A)/2

    lo_shift = (data[1, 0] - data[0, 0]) / 2
    la_shift = (data[3, 1] - data[0, 1]) / 2

    data[:, 0] -= lo_shift
    data[:, 1] -= la_shift
    s = 'POLYGON(('
    #     NOTE the double brakets, apparently necessary.

    for el in data:
        s += f'{el[0]} {el[1]},'

    # close the polygon by adding the first point again
    s += f'{data[0, 0]} {data[0, 1]}))'
    return s


def postgis_geom(data):
    point = f'POINT({data[0, 0]} {data[0, 1]})'

    poly = postgis_rect(data)
    return point, poly