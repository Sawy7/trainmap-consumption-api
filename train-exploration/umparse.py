import pandas as pd
import numpy as np
from datetime import datetime
import json
from tconsumption import calc_distance_two_points

def calculate_distance(row):
    return calc_distance_two_points((row["gps_longitude"], row["gps_latitude"]), (row["prev_gps_longitude"], row["prev_gps_latitude"]))

def um_csv_parser(csv_path):
    # Load from csv
    df = pd.read_csv(csv_path, delimiter=",")
    print(df.columns)

    # New cols (calculated)
    # Time
    df.insert(0, "datetime", pd.to_datetime(df["time"], unit="s"))
    del df["time"]
    df["delta_T"] = df["datetime"].diff(periods=1).abs()

    # Velocity
    df["velocity"] = np.sqrt(df["velocity_north"]**2+df["velocity_east"]**2)

    # Distance
    df["prev_gps_longitude"] = df["gps_longitude"].shift(1)
    df["prev_gps_latitude"] = df["gps_latitude"].shift(1)
    df["distance_m"] = df.apply(calculate_distance, axis=1)
    del df["prev_gps_longitude"]
    del df["prev_gps_latitude"]

    # Cumulative distance
    df['cumulative_distance_m'] = df['distance_m'].cumsum()

    return df

def make_geojson(lat, long):
    geojson = {
        "type": "LineString",
    }
    geojson["coordinates"] = list(zip(long, lat))
    return json.dumps(geojson)