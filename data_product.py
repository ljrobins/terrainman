
import src.terrainman as dp
import numpy as np
import matplotlib.pyplot as plt
import datetime
import os
import pyvista as pv
import pyspaceaware as ps

import numpy as np
import terrainman as tm
import matplotlib.pyplot as plt
station_lat_deg, station_lon_deg = 27.7172, 85.3240 # Katmandu
tdh = tm.TerrainDataHandler()
tile = tdh.load_tiles_containing(station_lat_deg, station_lon_deg)

sz, deg_radius = 3000, 1.0
lat_space = (station_lat_deg + deg_radius) - np.linspace(0, 2*deg_radius, sz)
lon_space = (station_lon_deg - deg_radius) + np.linspace(0, 2*deg_radius, sz)
lat_grid, lon_grid = np.meshgrid(lat_space, lon_space)
elev_grid = tile.interpolate(lat_grid, lon_grid) / 1e3

plt.imshow(elev_grid)
plt.show()

station_lat_deg, station_lon_deg = 27.7172, 85.3240 # Katmandu
# station_lat_deg, station_lon_deg = 43.65311150689344, -70.19252101245867 # Maine
station_lat_rad, station_lon_rad, station_alt_km = np.deg2rad(station_lat_deg), np.deg2rad(station_lon_deg), 0.00
station_height_above_terrain_km = 0.0

tdh = dp.TerrainDataHandler()
tile = tdh.load_tiles_containing(station_lat_deg, station_lon_deg)

sz = 3000
deg_radius = 1.0
lat_space = (station_lat_deg + deg_radius) - np.linspace(0, 2*deg_radius, sz)
lon_space = (station_lon_deg - deg_radius) + np.linspace(0, 2*deg_radius, sz)
lat_grid, lon_grid = np.meshgrid(lat_space, lon_space)
elev_grid = tile.interpolate(lat_grid, lon_grid) / 1e3

ps.tic()
station_ecef = ps.lla_to_itrf(station_lat_rad, station_lon_rad, station_alt_km)
ecef = ps.lla_to_itrf(np.deg2rad(lat_grid), np.deg2rad(lon_grid), elev_grid)
ll = np.hstack((lat_grid.reshape((lat_grid.size, 1)), lon_grid.reshape((lon_grid.size,1))))
station_terrain_alt_km = tile.interpolate(station_lat_deg, station_lon_deg) / 1e3
station_ecef = ps.lla_to_itrf(station_lat_rad, station_lon_rad, station_terrain_alt_km + station_height_above_terrain_km)
ps.toc()
print(station_terrain_alt_km)


station_enu = (ps.ecef_to_enu(station_ecef) @ station_ecef.T).T
enu = (ps.ecef_to_enu(station_ecef) @ ecef.T).T

x_grid = enu[:,0].reshape(lat_grid.shape) - station_enu[0,0]
y_grid = enu[:,1].reshape(lat_grid.shape) - station_enu[0,1]
z_grid = enu[:,2].reshape(lat_grid.shape) - station_enu[0,2]

elev_grid[elev_grid == 0] = np.nan
dem = pv.StructuredGrid(x_grid, y_grid, z_grid)
dem["Elevation [km]"] = elev_grid.flatten(order="F")
dem["Latitude"] = lat_grid.flatten(order="F")
dem["Longitude"] = lon_grid.flatten(order="F")

pl = pv.Plotter()
pl.add_mesh(dem, smooth_shading=True, scalars="Elevation [km]")
ps.scatter3(pl, np.array([[0,0,0]]), color='r', point_size=10)
# ps.plot_earth(pl, mode="ecef", stars=True)
pl.camera.focal_point = (0.0,0.0,0.0)
pl.show()

class AerosolDataHandler(dp.DataProduct):
    def __init__(self) -> None:
        # Signature (dates, cadence)
        self._storage_dir = os.path.join(os.environ['TERRAIN_DATA'], 'aerosol')
        self.url_base = 'https://ladsweb.modaps.eosdis.nasa.gov'
        self.download_format = "nc"
        self.save_format = "nc"
        self._VIIRS_AEROSOL_CADENCES = ["daily", "monthly"]

        self._after_init()
    
    def _set_fnames(self, date: datetime.datetime, cadence: str) -> str:
        assert cadence in self._VIIRS_AEROSOL_CADENCES, f"Cadence must be in {self._VIIRS_AEROSOL_CADENCES}!"
        csv_path = os.path.join(os.environ['TERRAIN_DATA'], f'viirs_aerosol_{cadence}.csv')
        with open(csv_path, 'r') as f:
            csv_els = f.read().split(',')
            if cadence == "daily":
                fnames = [c for c in csv_els if f'/{date.year}/{date.strftime("%j")}/' in c][-1]
            elif cadence == "monthly":
                month_first = datetime.datetime(year=date.year, month=date.month, day=1)
                month_first_doy = int(month_first.strftime("%j"))
                possible_date_strs = [f'/{date.year}/{d}/' for d in range(month_first_doy, month_first_doy+32)]
                fnames = [c for c in csv_els if any(map(c.__contains__, possible_date_strs))]
        self.extracted_fname = fnames.split('/')[-1]
        self.url_fname = fnames
        return self.extracted_fname

    def download(self, date: datetime.datetime, cadence: str):
        self._download_from_args((date, cadence))
    
    def load(self, date: datetime.datetime, cadence: str):
        return self._load_from_args((date, cadence))[0]


date = datetime.datetime(2020, 12, 9)   
cadence = "daily"
adh = AerosolDataHandler()
adh.download(date, cadence)
aset = adh.load(date, cadence)
adata = aset.variables['Aerosol_Optical_Thickness_550_Ocean_Maximum'][:]

# plt.imshow(adata)
# plt.show()

