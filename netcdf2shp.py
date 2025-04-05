import netCDF4 as nc
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import os
import argparse
from pyproj import CRS

def extract_netcdf_data(netcdf_file):
    """
    Extracts (lat, lon) and time-series data from a NetCDF file
    for variables that have dimensions (time, lat, lon).
    """
    dataset = nc.Dataset(netcdf_file, "r")

    # Check for coordinate variables
    if "lat" not in dataset.variables or "lon" not in dataset.variables:
        raise ValueError("NetCDF file must contain 'lat' and 'lon' variables.")

    latitudes = dataset.variables["lat"][:]
    longitudes = dataset.variables["lon"][:]

    # Filter variables with (time, lat, lon) dimensions
    valid_vars = [
        var for var in dataset.variables
        if hasattr(dataset.variables[var], "dimensions")
        and dataset.variables[var].ndim == 3
        and ("lat" in dataset.variables[var].dimensions and "lon" in dataset.variables[var].dimensions)
    ]

    if not valid_vars:
        raise ValueError("No 3D variables with (time, lat, lon) found in this NetCDF file.")

    records = []
    for i, lat in enumerate(latitudes):
        for j, lon in enumerate(longitudes):
            geom = Point(lon, lat)
            record = {"geometry": geom, "latitude": float(lat), "longitude": float(lon)}
            for var in valid_vars:
                try:
                    record[var] = dataset.variables[var][:, i, j].tolist()
                except Exception as e:
                    print(f"Warning: Skipped variable {var} due to error: {e}")
            records.append(record)

    dataset.close()
    return records

def save_as_shapefile(records, output_shp):
    """
    Saves extracted NetCDF data as a Shapefile using raw WKT to avoid EPSG lookup.
    
    Parameters:
        records (list): List of dictionaries containing extracted data.
        output_shp (str): Path to save the shapefile.
    """
    
    gdf = gpd.GeoDataFrame(records, geometry="geometry")

    # EPSG:4326 in raw WKT (hardcoded for full bypass)
    wkt_4326 = (
        'GEOGCRS["WGS 84",DATUM["World Geodetic System 1984",'
        'ELLIPSOID["WGS 84",6378137,298.257223563]],'
        'CS[ellipsoidal,2],AXIS["latitude",north],AXIS["longitude",east],'
        'UNIT["degree",0.0174532925199433]]'
    )

    gdf.set_crs(wkt_4326, inplace=True)  # Avoid EPSG lookup entirely
    os.makedirs(os.path.dirname(output_shp), exist_ok=True)
    gdf.to_file(output_shp)
    print("✓ Shapefile saved:", output_shp)

def save_as_csv(records, output_csv):
    df = pd.DataFrame(records)
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    df.to_csv(output_csv, index=False)
    print("✓ CSV saved:", output_csv)

def netcdf_to_geodata(netcdf_file, output_file, file_type="shp", verbose=False):
    records = extract_netcdf_data(netcdf_file)

    if file_type == "shp":
        save_as_shapefile(records, output_file)
    elif file_type == "csv":
        save_as_csv(records, output_file)
    else:
        raise ValueError("Invalid file_type. Choose 'shp' or 'csv'.")

    if verbose:
        print(f"Finished processing NetCDF: {netcdf_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert NetCDF to Shapefile or CSV.")
    parser.add_argument("--input", required=True, help="Path to the NetCDF file")
    parser.add_argument("--output", required=True, help="Path to the output file")
    parser.add_argument("--type", choices=["shp", "csv"], required=True, help="Output format")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

    args = parser.parse_args()
    netcdf_to_geodata(args.input, args.output, args.type, verbose=args.verbose)
