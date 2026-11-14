# Sealed handover
All supplied numbers are decimal exact values, not binary floating measurements.
The main archive covers 2026-10-20 and 2026-10-21. Hour 0 is the first day at 00:00.
The source logger exports retained only the rows in readings.csv; absent hours are missing, not zero.
Calibration sheet 17 certified the constant upper flow rates in meters.csv for this entire interval.
The display quantization intervals were already conservatively closed by metrology staff;
half_width is the full allowed distance from display, and includes both endpoints.
Logger clock alignment has been completed. There are no resets, wraps, leaks or storage tanks.
T is upstream of the three branches; these branches exhaust its outflow.
Valve-log 07 records the old routing; entry 08 records a completed instantaneous switch at hour 30.
The meter stays on its physical branch through the switch; the region receiving that branch changes.
All archive filenames and record IDs are stable source locators. origin preserves the upstream locator.

data/toy is a scaled two-period control archive, with period lengths of one hour.
data/conflict is another sealed control archive, differing only in the final trunk reading.
All control archive readings have zero half_width. No other hidden data is needed.
