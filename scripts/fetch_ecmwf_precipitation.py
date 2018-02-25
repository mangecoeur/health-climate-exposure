#!/usr/bin/env python
from pathlib import Path

import pandas as pd
from ecmwfapi import ECMWFDataServer

SETTINGS = {

}


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
        fetch_precipitation_year(out_folder, year)


if __name__ == '__main__':
    start_year = 1980
    end_year = 2018
    out_folder = Path('./monthly_precipitation')
    fetch_mean_precipitation(out_folder, start_year, end_year)
