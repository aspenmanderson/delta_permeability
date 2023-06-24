Process for inporting data from Excell sheet:

1) Add data
2) Right slick on the sheet, Display XY data (X: longitude, Y: latitude)
3) Export data as layer (should save features)
4) Data Management Tools > Projections and Transformations > Project
	- Input coridinate system: Geographic Cordinate System > World > WGS_1994 (appears as D_WGS_1994)
	- Out coridinate system: Projected Cordinate System> World > GCS_WGS_1994
5) Save projected data as layer (.shp file in Original Data Shapefiles)
	- If the projection works but the data does not appear in the correct spot, try porjecting again with the already projected data


(AA 6/24/2020)