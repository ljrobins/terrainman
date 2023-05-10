import src as pt
import numpy as np
import matplotlib.pyplot as plt
import os
import datetime
import io

from http.cookiejar import CookieJar
import urllib.request
import base64
from netCDF4 import Dataset
import pyspaceaware as ps

_VIIRS_AEROSOL_CADENCES = ["daily", "monthly"]
_VIIRS_URL_BASE = 'https://ladsweb.modaps.eosdis.nasa.gov'
def _viirs_files_from_date(date: datetime.datetime, cadence: str = "daily") -> list[str]:
    assert cadence in _VIIRS_AEROSOL_CADENCES, f"Cadence must be in {_VIIRS_AEROSOL_CADENCES}!"
    with open(f'src/data/viirs_aerosol_{cadence}.csv', 'r') as f:
        csv_els = f.read().split(',')
        if cadence == "daily":
            fnames = [c for c in csv_els if f'/{date.year}/{date.strftime("%j")}/' in c]
        elif cadence == "monthly":
            month_first = datetime.datetime(year=date.year, month=date.month, day=1)
            month_first_doy = int(month_first.strftime("%j"))
            possible_date_strs = [f'/{date.year}/{d}/' for d in range(month_first_doy, month_first_doy+32)]
            fnames = [c for c in csv_els if any(map(c.__contains__, possible_date_strs))]
        return fnames

def _viirs_download(fnames) -> str:
    extracted_fnames = [fname.split('/')[-1] for fname in fnames]
    for fname, extracted_fname in zip(fnames, extracted_fnames):
        url = f'{_VIIRS_URL_BASE}{fname}'
        cj = CookieJar()
        credentials = ('%s:%s' % (os.environ['EARTHDATA_USERNAME'], os.environ['EARTHDATA_PASSWORD']))
        encoded_credentials = base64.b64encode(credentials.encode('ascii'))
        req = urllib.request.Request(url, None, {'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8','Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3','Accept-Encoding': 'gzip, deflate, sdch','Accept-Language': 'en-US,en;q=0.8','Connection': 'keep-alive', 'User-Agent': u'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36',
                                                'Authorization': 'Basic %s' % encoded_credentials.decode("ascii")})

        opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cj))
        response = opener.open(req)
        bytesio_object = io.BytesIO(response.read())
        with open(extracted_fname, "wb") as f:
            f.write(bytesio_object.getbuffer())
    return extracted_fnames

def viirs_load(fname: str, variable_names: str | list[str], zero_fill: bool = True) -> list[np.ndarray]:
    variable_names = [variable_names] if variable_names is str else variable_names
    data_dict = {}
    rootgrp = Dataset(fname, "r", format="NETCDF4")
    for v in variable_names:
        data_dict[v] = rootgrp.variables[v][:].filled(0) if zero_fill else rootgrp.variables[v][:]
    rootgrp.close()
    return data_dict
        
if __name__ == "__main__":
    # lon, lat = 105, 43
    # lon = np.array([-83])
    # lat = np.array([35])
    # fpaths = pt.download_terrain_tiles_containing(lat, lon)
    # data = pt.load_terrain_tiles_containing(lat, lon)
    # elev = data[-1]
    # print(elev.shape, elev.dtype, elev.itemsize)
    # print(elev.nbytes / 1e6)

    # plt.imshow(elev)
    # plt.show()

    date = datetime.datetime(2020, 12, 9)
    fnames = _viirs_files_from_date(date)
    efnames = _viirs_download(fnames)
    x = viirs_load(efnames[0], variable_names=['Aerosol_Optical_Thickness_550_Land_Mean'])
    print(x)


#     print()
#     plt.imshow(val_land + val_ocean)
#     plt.show()

# 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/5200/AERDB_D3_VIIRS_SNPP/2023/098/AERDB_D3_VIIRS_SNPP.A{2023}098.002.{year}102000651.nc'
