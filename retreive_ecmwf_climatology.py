#!/usr/bin/env python

"""
Retieving climatology fields for ERA

Lancet countdown requirements 
-----------------------------

Heat exposure: Summer mean temperatures

Heatwaves:  99th percentile of the daily minimum temperature 



ERA climatology references
--------------------------

ERA Interim daily climatology (ERAI­DACL) consists of horizontal fields of statistical characteristics of 
atmospheric parameters at a given time and day of a year. For a given parameter and level and a given day of a year 
and time there are fields of 

- mean (a field of the filtered mean values of the parameter valid at a given hour and on a day of a year computed 
from the sampling 20­years period), 

- standard deviation (idem for standard deviation),

- minima and maxima (idem for extreme values)

- quantiles (the terciles and selected percentiles of the distribution of the parameter during the sampling period).

The ERAI­DACL fields are available from ECMWF MARS in the GRIB­1 format. The upper­air parameters are available at 00 
and 12 UTC network times and the surface and screen­level parameters at times 00, 06, 12, 18 UTC. Due to their high 
computational cost and data size the quantiles are available only for selected parameters. Climatology for wave 
parameters is available as well using different STREAM DACW 

The list of parameters and MARS identifiers

class=EI
stream=DACL (or DACW for wave parameters swh and mwp)
expver=0001
type=EM(mean)/ES(standard deviation)/FCMIN(minima)/FCMAX(maxima)
  for time=00:00/12:00
    levtype=PL 
        levelist=1/2/3/5/7/10/20/30/50/70/100/150/200/250/300/400/500/600/700/850/925/1000 
        param=Z/T/U/V/WIND/W/R/Q/O3
    levtype=PT 
        levelist=300/315/330/350/370/395 
        param=PV
    levtype=PV 
        levelist=2000
        param=PT
  for time=00:00/06:00/12:00/18:00
    levtype=SFC 
    param=2T/10U/10V/10SI/MSL/SSTK/CI/STL1/ISTL1/SWH/MWP
type=PB(quantiles) 
quantile=1:100/2:100/5:100/10:100/15:100/20:100/25:100/30:100/1:3/35:100/40:100/45:100/50:100/ 55:100/60:100/65:100/2:3/70:100/75:100/80:100/85:100/90:100/95:100/98:100/99:100
  At time 00:00/12:00:
    levtype=PL
      levelist=50/100/200/250/500/700/850/1000, param=Z/T
      levelist=200/250/500/700/850, param=U,V,Q
      levelist=850, param=WIND
    levtype=SFC
        param=2T/10U/10V/10SI/MSL/SWH/MWP
to specify date choose any leap year (e.g.1994

"""


from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

# Mean temperature
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
    "target": "output1.grib",
})
#
# # Minimum temperatures
# server.retrieve({
#     "dataset": "interim",
#
#     "class": "ei",
#     "stream": "dacl",
#     "expver": "1",
#     "type": "fcmin",
#     "levtype": "sfc",
#     "time": "00:00:00/06:00:00/12:00:00/18:00:00",
#     "param": "2T",
#     "date": 1996,
#     "target": "output2.grib",
# })
#
# # Maximum temperatures
# server.retrieve({
#     "dataset": "interim",
#
#     "class": "ei",
#     "stream": "dacl",
#     "expver": "1",
#     "type": "fcmax",
#     "levtype": "sfc",
#     "time": "00:00:00/06:00:00/12:00:00/18:00:00",
#     "param": "2T",
#     "date": 1996,
#
#     "target": "output3.grib",
# })
