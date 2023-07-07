WITH points AS (
	SELECT
		ogc_fid,
		ST_DumpPoints(wkb_geometry) AS point,
		ST_DumpPoints(ST_Transform(wkb_geometry, 3035)) AS pointdtm
	FROM tram_testing_data
)
SELECT ogc_fid,
ST_AsGeoJSON(ST_MakeLine(ST_Translate(
	ST_Force3DZ((point).geom),
	0::double precision,
	0::double precision,
	ST_Value(dtm.rast,(pointdtm).geom)
) ORDER BY ((point).path))) AS geom
FROM points
LEFT JOIN dtm_eu AS dtm ON ST_Intersects(dtm.rast, (pointdtm).geom)
WHERE ogc_fid = 1
GROUP BY ogc_fid;