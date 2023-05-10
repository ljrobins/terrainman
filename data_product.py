
import src as dp
import numpy as np
import matplotlib.pyplot as plt
import datetime

lat, lon = 43, -106
lon_arr = np.array([-83, 0])
lat_arr = np.array([35, 0.0])

tdh = dp.TerrainDataHandler()
tdh.download_tile(lat, lon)
tdh.download_tiles_containing(lat_arr, lon_arr)
tdata = tdh.load_tile(lat, lon)

plt.imshow(tdata)
plt.show()

date = datetime.datetime(2020, 12, 9)   
cadence = "daily"
adh = dp.AerosolDataHandler()
adh.download(date, cadence)
aset = adh.load(date, cadence)
adata = aset.variables['Aerosol_Optical_Thickness_550_Ocean_Maximum'][:]

plt.imshow(adata)
plt.show()

