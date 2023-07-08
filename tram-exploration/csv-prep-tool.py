import os, argparse
from subprocess import run, PIPE
import pandas as pd
import numpy as np
import json
import psycopg2

def remove_leader(path: str):
    return path.split(".")[-1]

def convert_sig_unsig(number):
    if number < 32767:
        return number
    else:
        return number-65535

def tram_csv_parser(csv_path):
    # Load from csv
    df = pd.read_csv(csv_path, delimiter=";", skiprows=73, decimal=",")

    # Rename cols
    rename_dict = {}
    for col in df.columns:
        new_name = remove_leader(col)
        rename_dict[col] = new_name
    df.rename(columns=rename_dict, inplace=True)

    # New cols (calculated)
    # Time
    df.insert(2, "DateTime", pd.to_datetime(df["Date"] + ' ' + df["Time"], format="%d.%m.%Y %H:%M:%S:%f"))
    del df["Date"]
    del df["Time"]
    df["delta_T"] = df["DateTime"].diff(periods=-1).abs()

    # Trakcni menic 1
    df["tm_stejnoproud_trakce"] = df["AM_A_C1_14_IwTCU_ILF1"].apply(convert_sig_unsig)/10
    df["tm_stejnoproud_spotrebice"] = df["AM_A_C1_14_IwAUX_I_DClink"]*0.1
    df["tm_stridaproud_vystupu"] = df["AM_A_C1_14_IwAUX_I_rms"]*0.01
    df["tm_vstupni_napeti"] = df["AM_A_C1_14_IwTCU_UD1Voltage"]/10
    df["tm_tah_vyvijeny_motory"] = df["AM_B_C1_15_IwTCU_E_MC1"].apply(convert_sig_unsig)/100
    df["tm_rychlost_3_napravy"] = df["VrTCU_A_A1_Velocity"]/3.6
    df["tm_prikon_trakce"] = df["tm_stejnoproud_trakce"]*df["tm_vstupni_napeti"]/1000
    df["tm_vykon_trakce"] = df["tm_tah_vyvijeny_motory"]*df["tm_rychlost_3_napravy"]
    df["tm_prikon_spotrebice"] = df["AM_A_C1_14_IwTCU_P_AUX"]/100
    df["tm_prikon_topeni_klima"] = 1.73*400*df["tm_stridaproud_vystupu"]/1000
    df["tm_vykon_zmareny"] = np.where((df["tm_vykon_trakce"] > 0) & (df["tm_vykon_trakce"] < df["tm_prikon_trakce"]), df["tm_prikon_trakce"]-df["tm_vykon_trakce"], 0)

    # Trakcni menic 2
    df["tm2_stejnoproud_trakce"] = df["AM_B_C1_15_IwTCU_ILF1"].apply(convert_sig_unsig)/10
    df["tm2_vstupni_napeti"] = df["AM_B_C1_15_IwTCU_UD1Voltage"]/10
    df["tm2_tah_vyvijeny_motory"] = df["AM_B_C1_15_IwTCU_E_MC1"].apply(convert_sig_unsig)/100
    df["tm2_rychlost_5_napravy"] = df["VrTCU_B_A1_Velocity"]/3.6
    df["tm2_prikon_trakce"] = df["tm2_stejnoproud_trakce"]*df["tm2_vstupni_napeti"]/1000
    df["tm2_vykon_trakce"] = df["tm2_tah_vyvijeny_motory"]*df["tm2_rychlost_5_napravy"]
    df["tm2_vykon_zmareny"] = np.where((df["tm2_vykon_trakce"] > 0) & (df["tm2_vykon_trakce"] < df["tm2_prikon_trakce"]), df["tm2_prikon_trakce"]-df["tm2_vykon_trakce"], 0)

    # Cela tramvaj
    df["cela_vykon_spotrebovany"] = np.where(df["tm_prikon_trakce"]+df["tm_prikon_spotrebice"]+df["tm2_prikon_trakce"] > 0, df["tm_prikon_trakce"]+df["tm_prikon_spotrebice"]+df["tm2_prikon_trakce"], 0)
    df["cela_vykon_spotrebovany_bez_spotrebicu"] = np.where(df["tm_prikon_trakce"]+df["tm2_prikon_trakce"] > 0, df["tm_prikon_trakce"]+df["tm2_prikon_trakce"], 0)
    df["cela_vykon_rekuperovany"] = np.where(df["tm_prikon_trakce"]+df["tm_prikon_spotrebice"]+df["tm2_prikon_trakce"] < 0, -(df["tm_prikon_trakce"]+df["tm_prikon_spotrebice"]+df["tm2_prikon_trakce"]), 0)
    df["cela_vykon_rekuperovany_bez_spotrebicu"] = np.where(df["tm_prikon_trakce"]+df["tm2_prikon_trakce"] < 0, -(df["tm_prikon_trakce"]+df["tm2_prikon_trakce"]), 0)
    df["cela_vykon_zmareny"] = df["tm_vykon_zmareny"]+df["tm2_vykon_zmareny"]
    df["cela_vykon_trakce"] = df["tm_vykon_trakce"]+df["tm2_vykon_trakce"]
    df["cela_prikon_topeni_klima"] = df["tm_prikon_topeni_klima"]

    df["cela_energie_spotrebovana"] = (df["cela_vykon_spotrebovany"]*df["delta_T"].dt.total_seconds()/3600).cumsum().shift(1)
    df.loc[0, "cela_energie_spotrebovana"] = 0.0

    df["cela_energie_spotrebovana_bez_spotrebicu"] = (df["cela_vykon_spotrebovany_bez_spotrebicu"]*df["delta_T"].dt.total_seconds()/3600).cumsum().shift(1)
    df.loc[0, "cela_energie_spotrebovana_bez_spotrebicu"] = 0.0

    df["cela_energie_rekuperovana"] = (df["cela_vykon_rekuperovany"]*df["delta_T"].dt.total_seconds()/3600).cumsum().shift(1)
    df.loc[0, "cela_energie_rekuperovana"] = 0.0

    df["cela_energie_rekuperovana_bez_spotrebicu"] = (df["cela_vykon_rekuperovany_bez_spotrebicu"]*df["delta_T"].dt.total_seconds()/3600).cumsum().shift(1)
    df.loc[0, "cela_energie_rekuperovana_bez_spotrebicu"] = 0.0

    df["cela_energie_zmarena"] = (df["cela_vykon_zmareny"]*df["delta_T"].dt.total_seconds()/3600).cumsum().shift(1)
    df.loc[0, "cela_energie_zmarena"] = 0.0

    df["cela_energie_topeni_klima"] = (df["cela_prikon_topeni_klima"]*df["delta_T"].dt.total_seconds()/3600).cumsum().shift(1)
    df.loc[0, "cela_energie_topeni_klima"] = 0.0

    df["km_total"] = df["QdwMesitTotalDistance"]*0.1

    return df

def make_geojson(lat, long):
    geojson = {
        "type": "LineString",
    }
    geojson["coordinates"] = list(zip(long, lat))
    return json.dumps(geojson)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", action="store", help=".csv file to parse", required=True)
    parser.add_argument("--host", action="store", help="DB server IP", required=True)
    parser.add_argument("--dbname", action="store", help="DB name", required=True)
    parser.add_argument("--dbuser", action="store", help="DB username", required=True)
    parser.add_argument("--dbpass", action="store", help="DB user password", required=True)
    parser.add_argument("--dbtable", action="store", help="DB table name", required=True)
    args = parser.parse_args()
    config = vars(args)

    df = tram_csv_parser(config["input"])
    tmp_file = "/tmp/tram-output.geojson"

    geojson = make_geojson(df["VrGpsLatitude"], df["VrGpsLongitude"])
    with open(tmp_file, "w") as f:
        f.write(geojson)

    p = run([
        'ogr2ogr', '-f', 'PostgreSQL',
        f'PG:host={config["host"]} dbname={config["dbname"]} user={config["dbuser"]} password={config["dbpass"]}',
        tmp_file, '-nln', config["dbtable"], '-append'
    ], stdout=PIPE)

    os.remove(tmp_file)
    
    conn = psycopg2.connect(f'dbname={config["dbname"]} user={config["dbuser"]} password={config["dbpass"]} host={config["host"]}')
    cur = conn.cursor()
    outputquery = f"""
    WITH points AS (
	    SELECT
	    	ogc_fid,
	    	ST_DumpPoints(wkb_geometry) AS point,
	    	ST_DumpPoints(ST_Transform(wkb_geometry, 3035)) AS pointdtm
	    FROM {config['dbtable']}
    )
    SELECT ST_AsGeoJSON(ST_MakeLine(ST_Translate(
	    ST_Force3DZ((point).geom),
	    0::double precision,
	    0::double precision,
	    ST_Value(dtm.rast,(pointdtm).geom)
    )    ORDER BY ((point).path))) AS geom
    FROM points
    LEFT JOIN dtm_eu AS dtm ON ST_Intersects(dtm.rast, (pointdtm).geom)
    WHERE ogc_fid = (SELECT MAX(ogc_fid) FROM {config['dbtable']})
    GROUP BY ogc_fid;
    """
    cur.execute(outputquery)
    rows = cur.fetchall()
    if len(rows) != 1:
        print("ERROR: Expected 1 row. Exiting...")

    dbjson = json.loads(rows[0][0])
    last_coord_idx = len(dbjson["coordinates"])-1
    modjson = {"type": dbjson["type"]}
    modjson["station_orders"] = [0, last_coord_idx]
    modjson["velocity_ways"] = [{"start":0, "end":last_coord_idx, "velocity":50}]
    modjson["coordinates"] = dbjson["coordinates"]

    new_file = os.path.splitext(config["input"])[0] + ".geojson"
    with open(new_file, "w") as f:
        json.dump(modjson, f)
    cur.close()
