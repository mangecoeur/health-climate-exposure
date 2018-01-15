
Notes on country specific data
==============================

Exposure values
---------------

IMPORTANT: DO NOT PERFORM AVERAGES OF EXPOSURE TIME SERIES ACROSS COUNTRIES

The mean exposure values represent the average of the given indicator 
(temperature change, length of heatwave, labour capacity) for a country
weighted by the population distribution within that country on a 
gridded basis.

So 

Sum(pop_i * indicator_i) / P

where pop_i is the population in a given grid cell, indicator_i is the value
of the indicator in that cell, P is the total population for the selected area.

Calculating the weighted average over the area of several selected countries 
will not therefore give you the same result as calculating the average of the
individual country results (because the normalisation by population total will
be wrong)



