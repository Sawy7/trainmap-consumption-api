import json
import psycopg2

def get_rail(config, relcislo, reversed):
    conn = psycopg2.connect(f'dbname={config["dbname"]} user={config["dbuser"]} password={config["dbpass"]} host={config["host"]}')
    cur = conn.cursor()
    outputquery = f"""
    SELECT ST_AsGeoJSON(ST_Collect(geom)) AS geojson FROM (
    SELECT (ST_DumpPoints(geom)).geom FROM even_processed_routes_line_dtm
    WHERE relcislo = {relcislo}
    ORDER BY (ST_DumpPoints(geom)).path[1]
    """
    if reversed:
        outputquery += " DESC"
    outputquery += ") AS all_points;"

    cur.execute(outputquery)
    rows = cur.fetchall()
    if len(rows) != 1:
        print("ERROR: Expected 1 row. Exiting...")

    dbjson = json.loads(rows[0][0])
    modjson = {"type": dbjson["type"]}
    modjson["coordinates"] = dbjson["coordinates"]

    # Stations
    outputquery = f"""
    SELECT even_station_relation.station_order
    FROM even_station_relation JOIN
    all_stations ON even_station_relation.station_id = all_stations.id
    WHERE relcislo = {relcislo}
    ORDER BY even_station_relation.relcislo, station_order
    """
    if reversed:
        outputquery += " DESC"

    cur.execute(outputquery)
    rows = cur.fetchall()
    modjson["station_orders"] = [x[0] for x in rows]
    if reversed:
        modjson["station_orders"] = [len(modjson["coordinates"])-1-x for x in modjson["station_orders"]]
    # Fix first and last station
    modjson["station_orders"][0] = 0
    modjson["station_orders"][-1] = len(modjson["coordinates"])-1

    # Velocity ways
    outputquery = f"""SELECT start_order, end_order, maxspeed
    FROM get_even_route_line_ways({relcislo}) AS ewr JOIN
    osm_ways ON ewr.way_id = osm_ways.id
    WHERE relcislo = {relcislo}
    """

    cur.execute(outputquery)
    rows = cur.fetchall()
    modjson["velocity_ways"] = [{"start":x[0], "end":x[1], "velocity":x[2]} for x in rows]

    new_file = config["output"]
    with open(new_file, "w") as f:
        json.dump(modjson, f)
    cur.close()
