import pandas as pd
import numpy as np
from datetime import datetime
import json
from tconsumption import calc_distance_two_points

def calculate_distance(row):
    return calc_distance_two_points((row["gps_longitude"], row["gps_latitude"]), (row["prev_gps_longitude"], row["prev_gps_latitude"]))

def calc_distance_all(df):
    df["prev_gps_longitude"] = df["gps_longitude"].shift(1)
    df["prev_gps_latitude"] = df["gps_latitude"].shift(1)
    df["distance_m"] = df.apply(calculate_distance, axis=1)
    del df["prev_gps_longitude"]
    del df["prev_gps_latitude"]

def fix_gps_outliers(df):
    max_distance = 1000
    to_fix = []
    cor_point = None

    for index, row in df.iterrows():
        if cor_point is not None:
            dist = calc_distance_two_points(
                (row["gps_longitude"], row["gps_latitude"]),
                (cor_point["gps_longitude"], cor_point["gps_latitude"])
            )
            if dist > max_distance:
                to_fix[-1]["points"].append(index)
            else:
                to_fix[-1]["cor_point2"] = row
                cor_point = None
        elif row["distance_m"] > max_distance:
            cor_point = df.iloc[index-1]
            to_fix.append({"cor_point": cor_point, "points": []})
            to_fix[-1]["points"].append(index)

    for batch in to_fix:
        lon_delta = batch["cor_point2"]["gps_longitude"] - batch["cor_point"]["gps_longitude"]
        lat_delta = batch["cor_point2"]["gps_latitude"] - batch["cor_point"]["gps_latitude"]
        gps_count = len([x for x in batch["points"] if df.iloc[x]["packet_type"] == "UM7GPSPacket"])
        lon_shift = 0
        lat_shift = 0
        for pidx in batch["points"]:
            current_point = df.iloc[pidx]
            if current_point["packet_type"] == "UM7GPSPacket":
                lon_shift += lon_delta/gps_count
                lat_shift += lat_delta/gps_count
            df.loc[pidx] = {
                "gps_longitude": batch["cor_point"]["gps_longitude"] + lon_shift,
                "gps_latitude": batch["cor_point"]["gps_latitude"] + lat_shift,
            }

def um_csv_parser(csv_path):
    # Load from csv
    df = pd.read_csv(csv_path, delimiter=",")
    print(df.columns)

    # New cols (calculated)
    # Time
    df.insert(0, "datetime", pd.to_datetime(df["time"], unit="s"))
    # del df["time"]
    df["delta_T"] = df["datetime"].diff(periods=1).abs()

    # Velocity
    # df["velocity"] = np.sqrt(df["velocity_north"]**2+df["velocity_east"]**2)

    # Distance
    calc_distance_all(df)

    # Fix GPS outliers
    fix_gps_outliers(df)

    # Recalc distance
    calc_distance_all(df)

    # Cumulative distance
    df['cumulative_distance_m'] = df['distance_m'].cumsum()

    return df

def make_geojson(lat, long):
    geojson = {
        "type": "LineString",
    }
    geojson["coordinates"] = list(zip(long, lat))
    return json.dumps(geojson)