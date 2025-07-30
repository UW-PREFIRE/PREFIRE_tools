"""
Define unit conversion parameters/factors
"""

# Larger to smaller units:
s_to_ms = 1.e3  # [s] -> [ms] (floating-point form)
int_s_to_ms = 1000  # [s] -> [ms] (integer form)

s_to_dms = 1.e4  # [s] -> [decimilliseconds] (floating-point form)

ms_to_us = 1.e3  # [ms] -> [us] (floating-point form)
int_ms_to_us = 1000  # [ms] -> [us] (integer form)

day_to_s = 86400.  # [days] -> [s] (floating-point form) Earth tropical day
hour_to_s = 3600.  # [hours] -> [s] (floating-point form)
minute_to_s = 60.  # [minutes] -> [s] (floating-point form)

km_to_m = 1.e3  # [km] -> [m]
m_to_mm = 1.e3  # [m] -> [mm]
m_to_um = 1.e6  # [m] -> [um]

# Smaller to larger units:
ms_to_s = 1.e-3  # [ms] -> [s] (floating-point form)

dms_to_ms = 1.e-1  # [decimilliseconds] to [ms]
dms_to_s = 1.e-4  # [decimilliseconds] to [s]

us_to_ms = 1.e-3  # [us] -> [ms] (floating-point form)

s_to_day = 1./86400.  # [s] -> [days] (floating-point form)
s_to_hour = 1./3600.  # [s] -> [hours] (floating-point form)
s_to_minute = 1./60.  # [s] -> [minutes] (floating-point form)

m_to_km = 1.e-3  # [m] -> [km]
mm_to_m = 1.e-3  # [mm] -> [m]
um_to_m = 1.e-6  # [um] -> [m]
