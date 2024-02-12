from tramparse import tram_csv_parser
from random import randint
import pandas as pd
import pyproj
import math
import json

DATA_PATH="../../enet-sz-data/real_rides/"

def calc_distance_two_points(point_a, point_b):
    R = 6371e3  # meters
    φ1 = point_a[1] * math.pi/180  # φ, λ in radians
    φ2 = point_b[1] * math.pi/180
    Δφ = (point_b[1]-point_a[1]) * math.pi/180
    Δλ = (point_b[0]-point_a[0]) * math.pi/180

    a = math.sin(Δφ/2) * math.sin(Δφ/2)  \
        + math.cos(φ1) * math.cos(φ2) * math.sin(Δλ/2) * math.sin(Δλ/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))

    return R * c  # meters

def calc_coordinates_from_dist(distance, point_a):
    R = 6371e3  # meters
    φ1 = point_a[1] * math.pi/180  # φ, λ in radians
    λ1 = point_a[0] * math.pi/180

    Δ = distance / R  # angular distance

    φ2 = math.asin(math.sin(φ1) * math.cos(Δ) + math.cos(φ1) * math.sin(Δ) * math.cos(0))
    Δλ = math.atan2(math.sin(0) * math.sin(Δ) * math.cos(φ1), math.cos(Δ) - math.sin(φ1) * math.sin(φ2))
    λ2 = λ1 + Δλ

    φ2 = φ2 * 180 / math.pi
    λ2 = λ2 * 180 / math.pi

    return λ2, φ2

def make_anon_csv(input_csv, input_geojson, output_csv):
    df = tram_csv_parser(input_csv)

    df_out = df[[
        "km_total",
        "VrGpsLatitude", "VrGpsLongitude",
        "tm_rychlost_3_napravy",
        "cela_vykon_spotrebovany_bez_spotrebicu", "cela_vykon_rekuperovany_bez_spotrebicu",
        "cela_energie_spotrebovana_bez_spotrebicu", "cela_energie_rekuperovana_bez_spotrebicu",
        "cela_vykon_vyrovnany_bez_spotrebicu", "cela_energie_vyrovnana_bez_spotrebicu"
    ]]
    df_out = df_out.rename(columns={
        "km_total": "celkem_km",
        "VrGpsLatitude": "gps_lat",
        "VrGpsLongitude": "gps_long",
        "tm_rychlost_3_napravy": "rychlost_ms"
    })

    distances = []
    for i in range(1, len(df)):
        lat1, lon1 = df_out.iloc[i-1]['gps_lat'], df_out.iloc[i-1]['gps_long']
        lat2, lon2 = df_out.iloc[i]['gps_lat'], df_out.iloc[i]['gps_long']
        distance = calc_distance_two_points((lon1, lat1), (lon2, lat2))
        distances.append(distance)

    # offset_lat = randint(1, 2)
    # offset_long = randint(1, 2)
    # prev_p = (df_out["gps_long"][0], df_out["gps_lat"][0])
    # df_out["gps_long"].iloc[0] = df_out["gps_long"].iloc[0] + offset_long
    # df_out["gps_lat"].iloc[0] = df_out["gps_lat"].iloc[0] + offset_lat
    # for i in range(1, len(df_out)):
    #     p = (df_out.iloc[i]["gps_long"], df_out.iloc[i]["gps_lat"])
    #     distance = calc_distance_two_points(prev_p, p)
    #     prev_p = (df_out["gps_long"][i], df_out["gps_lat"][i])
    #     if distance == 0:
    #         df_out["gps_long"].iloc[i] = df_out["gps_long"].iloc[i-1]
    #         df_out["gps_lat"].iloc[i] = df_out["gps_lat"].iloc[i-1]
    #     else:
    #         new_p = calc_coordinates_from_dist(
    #             distance,
    #             (df_out["gps_long"][i-1], df_out["gps_lat"][i-1])
    #         )
    #         df_out["gps_long"].iloc[i] = new_p[0]
    #         df_out["gps_lat"].iloc[i] = new_p[1]

    wgs84_crs = pyproj.CRS.from_epsg(4326)
    local_crs = pyproj.CRS.from_epsg(5514) 
    to_5514 = pyproj.Transformer.from_crs(wgs84_crs, local_crs, always_xy=True)
    to_gps = pyproj.Transformer.from_crs(local_crs, wgs84_crs, always_xy=True)

    df_out[["loc_x", "loc_y"]] = df_out.apply(
        lambda row: pd.Series(
            to_5514.transform(row["gps_long"], row["gps_lat"])
        ), axis=1)

    offset_x = randint(500, 1000)*1000
    offset_y = randint(500, 1000)*1000

    df_out["loc_x"] = df_out["loc_x"] + offset_x
    df_out["loc_y"] = df_out["loc_y"] + offset_y

    df_out[["gps_long", "gps_lat"]] = df_out.apply(
        lambda row: pd.Series(
            to_gps.transform(row["loc_x"], row["loc_y"])
        ), axis=1)

    distances2 = []
    for i in range(1, len(df)):
        lat1, lon1 = df_out.iloc[i-1]['gps_lat'], df_out.iloc[i-1]['gps_long']
        lat2, lon2 = df_out.iloc[i]['gps_lat'], df_out.iloc[i]['gps_long']
        distance = calc_distance_two_points((lon1, lat1), (lon2, lat2))
        distances2.append(distance)

    df_out = df_out.drop(["loc_x", "loc_y"], axis=1)

    # for i in range(len(distances)):
    #     print(i, distances[i], distances2[i])

    with open(input_geojson) as f:
        lines = f.readlines()
        geojson_raw = "".join(lines)
        geojson = json.loads(geojson_raw)
        elevation = [p[2] for p in geojson["coordinates"]]
        df_out.insert(loc=3, column="vyskovy_profil", value=elevation)

    df_out.to_csv(output_csv, index=False)

if __name__ == "__main__":
    make_anon_csv(
        DATA_PATH + "/DPO/Jízda_Poruba_Zátiší_20.04.2022/01-Vozovna-Vresinska_2022-04-20.csv",
        DATA_PATH + "/DPO/Jízda_Poruba_Zátiší_20.04.2022/01-Vozovna-Vresinska_2022-04-20.geojson",
        DATA_PATH + "/DPO-Anonymous/01.csv"
    )

    make_anon_csv(
        DATA_PATH + "DPO/Jízda_Poruba_Zátiší_20.04.2022/02-Vresinska-Zatisi_2022-04-20.csv",
        DATA_PATH + "DPO/Jízda_Poruba_Zátiší_20.04.2022/02-Vresinska-Zatisi_2022-04-20.geojson",
        DATA_PATH + "/DPO-Anonymous/02.csv"
    )

    make_anon_csv(
        DATA_PATH + "DPO/Jízda_Poruba_Zátiší_20.04.2022/03-Zatisi-Vresinska_2022-04-20.csv",
        DATA_PATH + "DPO/Jízda_Poruba_Zátiší_20.04.2022/03-Zatisi-Vresinska_2022-04-20.geojson",
        DATA_PATH + "/DPO-Anonymous/03.csv"
    )

    make_anon_csv(
        DATA_PATH + "DPO/Jizdy_Centrum_07-08-12_07.2022/1710_02.csv",
        DATA_PATH + "DPO/Jizdy_Centrum_07-08-12_07.2022/1710_02.geojson",
        DATA_PATH + "/DPO-Anonymous/04.csv"
    )