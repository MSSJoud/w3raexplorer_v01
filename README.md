# W3RA Explorer QGIS Plugin

**Author:** Mehdi S. Joud  
**Email:** mshjo@plan.aau.dk  
**Version:** .01  
**License:** MIT

---

## What It Does

The **W3RA Explorer** plugin helps you explore gridded hydrological data from the **W3RA model** stored in **NetCDF format**.

- Load NetCDF datasets of water storage variables (e.g. `Sg_EU_2010`, `Sd_EU_2011`, etc.)
- Convert them to **shapefile point layers** in QGIS
- Click on any point to **interactively view time series** and **anomaly plots**
- Apply **regression fits**: linear, polynomial (dynamic degree), Gaussian smoothing, Fourier
- Identify variables using **dropdowns**, not raw input
- **Export selected points to CSV** (only in the *complete version*) contact the author directly to request access.

---

## 🛠 Requirements

QGIS ≥ 3.16  
Python 3.10+  
Python modules used:
- numpy
- matplotlib
- netCDF4
- PyQt5
- geopandas
- shapely
- xarray
- scikit-learn (for regressions)
- pandas

> ⚠ Make sure the required packages are installed in the QGIS Python environment or use system-level packages via `apt`.

---

## Optional Preprocessing

If you have `.mat` files from the W3RA model, convert them to NetCDF with the script:

## Installation

1. Ensure your system has the following dependencies installed:

2. Open QGIS → Plugins → Manage and Install Plugins → Install from ZIP → Select `W3RAExplorer.zip`.

3. You’ll find the plugin button in the QGIS toolbar: **W3RA Explorer**

## Tools

This plugin works with NetCDF-formatted W3RA output.
To convert W3RA `.mat` files to NetCDF, use the script provided: tools/W3RA_mat_2netcdf.py
python3 tools/W3RA_mat_2netcdf.py --input_dir ... --output_file ...
