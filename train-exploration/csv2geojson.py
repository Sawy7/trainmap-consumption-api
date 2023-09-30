import os, argparse
from subprocess import run, PIPE
import pandas as pd
import numpy as np
import json
import psycopg2
from umparse import um_csv_parser, make_geojson

def prep(config, df=None):
    if df is None:
        df = um_csv_parser(config["input"])
    tmp_file = "/tmp/tram-output.geojson"

    geojson = make_geojson(df["gps_latitude"], df["gps_longitude"])
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
    modjson["coordinates"] = dbjson["coordinates"]

    # Velocity ways
    outputquery = f"""
    SELECT min(index) AS start_order, max(index) AS end_order, maxspeed FROM (
		SELECT *, count(is_reset) OVER (ORDER BY index) AS grp FROM (
			SELECT index, dpoints.geom, closest.id AS way_id, closest.maxspeed,
			CASE WHEN LAG(closest.id) OVER (ORDER BY index) <> closest.id THEN 1 END AS is_reset
			FROM
			(
				SELECT (ST_DumpPoints(wkb_geometry)).path[1]-1 AS index, ST_Force2D((ST_DumpPoints(wkb_geometry)).geom) AS geom
				FROM {config['dbtable']}
				WHERE ogc_fid = (SELECT MAX(ogc_fid) FROM {config['dbtable']})
			) AS dpoints
			JOIN LATERAL
			(
				SELECT id, maxspeed
				FROM osm_ways
				WHERE maxspeed IS NOT NULL
				ORDER BY dpoints.geom <-> osm_ways.geom
				LIMIT 1
			) AS closest
			ON true
		) AS indexed
	) AS grouped
	GROUP BY way_id, grp, maxspeed
	ORDER BY min(index);
    """
    cur.execute(outputquery)
    rows = cur.fetchall()
    modjson["velocity_ways"] = [{"start":x[0], "end":x[1], "velocity":x[2]} for x in rows]

    new_file = config["output"]
    with open(new_file, "w") as f:
        json.dump(modjson, f)
    cur.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", action="store", help=".csv file to parse", required=True)
    parser.add_argument("--output", action="store", help=".geojson file to save", required=True)
    parser.add_argument("--host", action="store", help="DB server IP", required=True)
    parser.add_argument("--dbname", action="store", help="DB name", required=True)
    parser.add_argument("--dbuser", action="store", help="DB username", required=True)
    parser.add_argument("--dbpass", action="store", help="DB user password", required=True)
    parser.add_argument("--dbtable", action="store", help="DB table name", required=True)
    args = parser.parse_args()
    config = vars(args)

    prep(config)
