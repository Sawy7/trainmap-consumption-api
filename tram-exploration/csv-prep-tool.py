import os, argparse
from subprocess import run, PIPE
import pandas as pd
import numpy as np
import json
import psycopg2
from tramparse import tram_csv_parser, make_geojson

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
