import scipy.io
import netCDF4 as nc
import numpy as np
import os
import subprocess


"""
W3RA_mat_2netcdf.py

Author: Mehdi Joud  
Created: 2025-04-04

Description:
------------
This script converts W3RA .mat output files into a consolidated NetCDF (.nc) file
that can be used as input to the W3RAExplorer QGIS plugin.

The .mat files contain yearly outputs from the W3RA hydrological model, including:
- Groundwater (Sg)
- Soil moisture layers (Ss, Sd, S0)
- Surface water (Sr)
- Snow water equivalent (Ssnow)

Features:
---------
- Extracts time-series data for each hydrological variable.
- Preserves variable names with yearly suffixes (e.g., Sg_EU_2010, Sg_EU_2011, ...).
- Consolidates multi-year datasets into a single NetCDF file.
- Compatible with the QGIS W3RAExplorer plugin for interactive point-based visualization.

Usage:
------
Place this script in your working directory with access to .mat files.
Run the script with:
    python3 W3RA_mat_2netcdf.py --input_dir /path/to/matfiles --output /path/to/output.nc

Dependencies:
-------------
- scipy
- xarray
- numpy
- os, argparse
"""

def convert_w3ra_mat_to_netcdf(input_dir, latlon_file, output_file):
    """
    Convert multiple W3RA .mat files into a single NetCDF file with separate layers for each year.

    Parameters:
        input_dir (str): Directory containing the .mat files.
        latlon_file (str): Path to the LatLon.mat file containing latitude and longitude grids.
        output_file (str): Path to save the combined NetCDF file.
    """
    # Load latitude and longitude from LatLon.mat
    latlon_data = scipy.io.loadmat(latlon_file, struct_as_record=False, squeeze_me=True)
    lat = np.array(latlon_data.get('lat_EU')).squeeze()
    lon = np.array(latlon_data.get('lon_EU')).squeeze()

    if lat is None or lon is None:
        raise ValueError(f"Latitude or Longitude variable missing in {latlon_file}. Available keys: {list(latlon_data.keys())}")

    # Get all .mat files and sort them by year
    mat_files = sorted([f for f in os.listdir(input_dir) if f.endswith(".mat") and "LatLon" not in f])
    
    # Extract years from filenames
    years = [int(f.split("_")[-1].split(".")[0]) for f in mat_files]  
    
    # Create NetCDF file
    with nc.Dataset(output_file, "w", format="NETCDF4") as ds:
        # Create dimensions
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))

        # Create coordinate variables
        lat_var = ds.createVariable("lat", "f8", ("lat",))
        lon_var = ds.createVariable("lon", "f8", ("lon",))

        # Assign CF-compliant attributes
        lat_var.units = "degrees_north"
        lon_var.units = "degrees_east"

        # Assign coordinate values
        lat_var[:] = lat
        lon_var[:] = lon

        # Process each file and store each year's data in separate variables
        for idx, mat_file in enumerate(mat_files):
            year = years[idx]
            print(f"Processing {year}: {mat_file}")

            # Load .mat file
            mat_data = scipy.io.loadmat(os.path.join(input_dir, mat_file), struct_as_record=False, squeeze_me=True)
            print(f"Variables in {mat_file}: {list(mat_data.keys())}")

            # Find the first valid data variable to determine days count
            data_keys = [key for key in mat_data.keys() if key.startswith('S') and isinstance(mat_data[key], np.ndarray)]
            if not data_keys:
                print(f"⚠️ No valid data variables found in {mat_file}, skipping year {year}")
                continue

            days_in_year = mat_data[data_keys[0]].shape[0]  # Extract from first valid variable
            
            # Create time dimension for each year
            ds.createDimension(f"time_{year}", days_in_year)
            time_var = ds.createVariable(f"time_{year}", "i4", (f"time_{year}",))
            time_var.units = f"days since {year}-01-01"
            time_var.calendar = "standard"
            time_var[:] = np.arange(1, days_in_year + 1)
            
            # Process each variable
            for var in data_keys:
                data = mat_data[var]
                print(f"Checking {var} - Available shape: {data.shape if isinstance(data, np.ndarray) else 'Not an array'}")

                if isinstance(data, np.ndarray) and data.ndim == 3:
                    var_name = f"{var}_{year}"  # Store each year as a separate layer
                    print(f"\u2705 Writing variable {var_name} to NetCDF with shape {data.shape}")

                    # Create variable for this year
                    data_var = ds.createVariable(var_name, "f4", (f"time_{year}", "lat", "lon"), zlib=True)
                    data_var.units = "mm"
                    data_var.long_name = f"{var} for {year}"
                    
                    # Assign data
                    data_var[:, :, :] = data
                else:
                    print(f"\u26A0 Skipping {var}_{year} due to shape mismatch or invalid data")

    print(f"Multi-year NetCDF saved at {output_file}")

# Example usage
if __name__ == "__main__":
    input_dir = "/home/ubuntu/data_AOI-3/w3ra_aoi_3/"
    latlon_file = "/home/ubuntu/data_AOI-3/w3ra_aoi_3/LatLon.mat"
    output_file = "/home/ubuntu/data_AOI-3/w3ra_aoi_3/netcdf/W3RA_2010_2024_.nc"

    convert_w3ra_mat_to_netcdf(input_dir, latlon_file, output_file)


# Extract all Sg_EU variables into a separate file 
    sg_output_file = output_file.replace("W3RA_2010_2024_.nc", "W3RA_Sg_EU.nc")
    extract_command = f"ncks -v $(ncdump -h {output_file} | grep Sg_EU | awk '{{print $1}}' | grep -o 'Sg_EU_[0-9]*' | tr '\n' ',') -O {output_file} {sg_output_file}"
    print(f"Running extraction command: {extract_command}")
    subprocess.run(extract_command, shell=True, check=True)
    print(f"\u2705 Extracted Sg_EU variables into {sg_output_file}")