import terrainman as tm

# terrainman needs to be fixed to query from, e.g., here:
# https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/SRTMGL1.003/N00E013.SRTMGL1.hgt/N00E013.SRTMGL1.hgt.zip

# as the old usgs data is now offline

tile = tm.TerrainDataHandler().load_tiles_containing(43.0, -110.0)