# Calculatate the wetbulb and  wetbulb globe temperature at each time step

# Tref = globe.tas[,, j]  # near surface air temperature
# rhref = globe.hurs[,, j]
# p = globe.psl[,, j] / 100
# esat = exp(
#     (-2991.2729 / Tref ^ 2) - (6017.0128 / Tref) + 18.87643854 - 0.028354721 * Tref + (1.7838301E-5) * Tref ^ 2 - (
#     8.4150417E-10) * Tref ^ 3 + (4.4412543E-13) * Tref ^ 4 + 2.858487 * log(Tref)) / 100
# wsat = 621.97 * esat / (p - esat)
# w = rhref / 100 * wsat
# TL = 1 / (1 / (Tref - 55) - log(rhref / 100) / 2840) + 55
# TE = Tref * (1000 / p) ^ (0.2854 * (1 - 0.28E-3 * w)) * exp((3.376 / TL - 0.00254) * w * (1 + 0.81E-3 * w))
# WBT[,, j]= WBT[,, j]+(45.114 - 51.489 * (TE / 273.15) ^ (-3.504))
# WBGT[,, j]= WBGT[,, j]+0.7 * (45.114 - 51.489 * (TE / 273.15) ^ (-3.504)) + 0.3 * (Tref - 273.15)

#  WBGT is the array of the sum of all model wet bulb globe temperatures at each time step

import numpy as np


# Tref = globe.tas[,, j]  # near surface air temperature
# rhref = globe.hurs[,, j]
# p = globe.psl[,, j] / 100
#


# TODO check whether this should be done in ˚C or Kelvin

def calculate_wbt(t_ref, relative_humidity, surface_pressure):
    """Empirical calculation of wet bulb temperature from temperature, humidity, and pressure

    Args:
        t_ref: Dry bulb near-surface air temperature (K)
        relative_humidity: Relative humidity (%)
        surface_pressure: Surface air pressure (Pa)

    Returns:
        Wet Bulb Temperature
        TODO: is this in ˚C or K??
    """

    # Empirical formula for e_sat
    a = (-2991.2729 / t_ref ** 2)
    b = -6017.0128 / t_ref
    c = -0.028354721 * t_ref
    d = 1.7838301E-5 * t_ref ** 2
    e = -8.4150417E-10 * t_ref ** 3
    f = 4.4412543E-13 * t_ref ** 4
    g = 2.858487 * np.log(t_ref)
    sat = 18.87643854 + a + b + c + d + e + f + g
    e_sat = np.exp(sat) / 100

    w_sat = 621.97 * e_sat / (surface_pressure - e_sat)
    humidity_frac = relative_humidity / 100
    w = humidity_frac * w_sat
    t_l = 1 / (1 / (t_ref - 55) - np.log(humidity_frac) / 2840) + 55
    t_e = t_ref * (1000 / surface_pressure) ** (0.2854 * (1 - 0.28E-3 * w)) * np.exp(
        (3.376 / t_l - 0.00254) * w * (1 + 0.81E-3 * w))

    wbt = 45.114 - 51.489 * (t_e / 273.15) ** (-3.504)

    return wbt


def calculate_wbgt(t_ref, relative_humidity, surface_pressure):
    """
    
    Args:
        t_ref: 
        relative_humidity: 
        surface_pressure: 

    Returns:
        WBGT (˚C)
    """
    wbt = calculate_wbt(t_ref, relative_humidity, surface_pressure)
    wbgt = 0.7 * wbt + 0.3 * (t_ref - 273.15)
    return wbgt
