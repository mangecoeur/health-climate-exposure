#!/usr/bin/env python
import argparse
from pathlib import Path

import pandas as pd
from ecmwfapi import ECMWFDataServer


def fetch_monthly_mean_temperature(out_folder, start=1980, end=2018):
    # TODO should retrive all months, not that big...
    # TODO convert to get a file per year
    server = ECMWFDataServer()
    years = range(start, end)
    for year in years:
        out_file = str(out_folder / "{year}_monthly_mean.nc".format(year=year))
        date_list = '/'.join((t.strftime('%Y%m%d') for t in pd.date_range('{year}-01-01'.format(year=year),
                                                                          '{year}-12-01'.format(year=year), freq='MS')))
        server.retrieve({
            "class": "ei",
            "dataset": "interim",
            "date": date_list,
            "expver": "1",
            "grid": "0.75/0.75",
            "levtype": "sfc",
            "param": "134.128/167.128/168.128",
            "stream": "moda",
            "type": "an",
            "format": 'netcdf',
            "target": out_file,
        })


def fetch_climatology(out_file):
    # Mean temperature
    server = ECMWFDataServer()

    server.retrieve({
        "dataset": "interim",
        "class": "ei",
        "stream": "dacl",
        "expver": "1",
        "type": "em",
        "levtype": "sfc",
        "time": "00:00:00/06:00:00/12:00:00/18:00:00",
        "param": "2T",
        "date": "1989-01-01/to/2005-01-01",
        "format": 'netcdf',
        "target": out_file,
    })



def fetch_single_summer_year(year, out_folder):
    server = ECMWFDataServer()
    # Select only DJF/JJA months
    rng = pd.date_range('{year}-01-01'.format(year=year),
                        '{year}-03-1'.format(year=year), freq='D').append(
        pd.date_range('{year}-06-01'.format(year=year),
                      '{year}-09-1'.format(year=year), freq='D')).append(
        pd.date_range('{year}-12-01'.format(year=year),
                      '{year}-12-31'.format(year=year), freq='D'))
    date_list = '/'.join((t.strftime('%Y-%m-%d') for t in rng))
    out_file = str(out_folder / "{year}_summer_temperature_2m.nc".format(year=year))
    server.retrieve({
        "class": "ei",
        "dataset": "interim",
        "date": date_list,
        "expver": "1",
        "grid": "0.75/0.75",
        "levtype": "sfc",
        "param": "167.128",
        "step": "0",
        "stream": "oper",
        "time": "00:00:00/06:00:00/12:00:00/18:00:00",
        "type": "an",
        "format": 'netcdf',
        "target": out_file,
    })
    return


def fetch_summer_temperatures(out_folder, start=1980, end=2018):
    """
    Get the 6-hourly temperatures

    Args:
        start:
        end:

    Returns:

    """

    for year in range(start, end):
        fetch_single_summer_year(year, out_folder)


def fetch_precipitation_year(out_folder, year):
    """
    Fetch the monthly mean of daily accumulation
    for total precipitation. Think of it as mean daily precipitation
    for each month (to get total PPT for month, need to multiply by length of month)

    Args:
        out_folder:
        year:

    Returns:

    """
    server = ECMWFDataServer()

    rng = pd.date_range('{year}-01-01'.format(year=year),
                        '{year}-12-31'.format(year=year), freq='M')
    date_list = '/'.join((t.strftime('%Y-%m-%d') for t in rng))
    out_file = str(out_folder / "{year}_precipitation.nc".format(year=year))

    server.retrieve({
        "class": "ei",
        "dataset": "interim",
        "date": date_list,
        "expver": "1",
        "grid": "0.75/0.75",
        "levtype": "sfc",
        "param": "228.128",
        "step": "0-12",
        "stream": "mdfa",
        "type": "fc",
        "target": out_file,
    })


def fetch_mean_precipitation(out_folder, start=1980, end=2018):
    for year in range(start, end):
        fetch_precipitation_year(year, out_folder)


# TODO CLI interface
def main():
    parser = argparse.ArgumentParser(description='Download ECMWF data.')
    parser.add_argument('out_folder', type=str)
    parser.add_argument('--variable', choices=['summer_temperature', 'monthly_precipitation'])
    parser.add_argument('--start', default=1980, type=int)
    # TODO change default to current year
    parser.add_argument('--end', default=2018, type=int)


if __name__ == '__main__':
    main()
    # summer_t_out = Path('~/Data/').expanduser() / 'weather' / 'ecmwf' / 'summer_temperature'
    # summer_t_out = Path('./summer_temperature')
    #
    # fetch_summer_temperatures(summer_t_out)

    # mean_t_out = Path('~/Data/').expanduser() / 'weather' / 'ecmwf' / 'monthly_means'
    # fetch_monthly_mean_temperature(mean_t_out)