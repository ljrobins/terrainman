import os
import pickle
import ssl
import urllib.request
import zipfile
import io
from abc import ABC
from typing import Union
import base64
from http.cookiejar import CookieJar
import numpy as np
import rasterio
from netCDF4 import Dataset
import shutil

def purge_data():
    assert 'terrainman' in os.environ["TERRAIN_DATA"]
    shutil.rmtree(os.environ["TERRAIN_DATA"])
    os.mkdir(os.environ["TERRAIN_DATA"])

def _strip_tuple(t: tuple) -> tuple:
    """Returns a copy of a tuple-of-tuples with no extra layers such that `t[0][0]` is not a tuple

    :param t: Input tuple
    :type t: tuple
    :return: Same tuple, stripped of outer layers past the first
    :rtype: tuple
    """
    while isinstance(t[0], tuple):
        t = t[0]
    return (t,)  # I am braindead and tired so this terribleness will have to do

class DataProduct(ABC):
    HTTP_ERR_MSG = "HTTP error occured, aborting..."
    DNE_MSG = "Data does not exist to download, skipping"
    EXISTS_MSG = "Data already downloaded, skipping"

    def _after_init(self):
        self._bad_input_path = os.path.join(self._storage_dir, "bad_inputs.pkl")
        if not os.path.exists(self._storage_dir):
            os.mkdir(self._storage_dir)
        if not os.path.exists(self._bad_input_path):
            with open(os.path.join(self._bad_input_path), "wb") as f:
                pickle.dump([], f)

    def _before_request(self, *args):
        try:
            os.environ["EARTHDATA_USERNAME"]
            os.environ["EARTHDATA_PASSWORD"]
        except KeyError as e:
            print(
                "\nEnvironmental variables 'EARTHDATA_USERNAME' and 'EARTHDATA_PASSWORD' must be set!\nRegister at https://urs.earthdata.nasa.gov/users/new and set them with either:\n    - os.environ['EARTHDATA_USERNAME'] = '<your_username>'\n      os.environ['EARTHDATA_PASSWORD'] = '<your_password>'\n    - Create a .env.shared file in the calling directory containing\nEARTHDATA_USERNAME=<your_username>\nEARTHDATA_PASSWORD=<your_password>"
            )
            raise
        efname = self._set_fnames(*args)
        if self._is_fname_bad(efname):
            print(self.DNE_MSG)
            return False
        if efname in os.listdir(self._storage_dir):
            return False
        print("File not found in storage, proceeding to download...")
        return True

    def _make_request(self, *args):
        print(f"Downloading {self.extracted_fname}...")
        self.request_url = f"{self.url_base}{self.url_fname}"
        print(f"Using: {self.request_url}")

        cj = CookieJar()
        credentials = "%s:%s" % (
            os.environ["EARTHDATA_USERNAME"],
            os.environ["EARTHDATA_PASSWORD"],
        )
        encoded_credentials = base64.b64encode(credentials.encode("ascii"))

        context = ssl.create_default_context()
        context.check_hostname = False
        context.verify_mode = ssl.CERT_NONE

        req = urllib.request.Request(
            self.request_url,
            None,
            {
                "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
                "Accept-Charset": "ISO-8859-1,utf-8;q=0.7,*;q=0.3",
                "Accept-Encoding": "gzip, deflate, sdch",
                "Accept-Language": "en-US,en;q=0.8",
                "Connection": "keep-alive",
                "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36",
                "Authorization": "Basic %s" % encoded_credentials.decode("ascii"),
            },
        )

        opener = urllib.request.build_opener(
            urllib.request.HTTPCookieProcessor(cj),
            urllib.request.HTTPSHandler(context=context),
        )

        try:
            self.response = opener.open(req)
        except urllib.error.HTTPError as e:
            print(self.HTTP_ERR_MSG)
            self._store_bad_input(*args)
            self.response = None

    def _after_request(self):
        if self.response is not None:
            if self.download_format == "zip":
                z = zipfile.ZipFile(io.BytesIO(self.response.read()))
                z.extractall(self._storage_dir)
                print("Done.")
            else:
                bytesio_object = io.BytesIO(self.response.read())
                with open(
                    os.path.join(self._storage_dir, self.extracted_fname), "wb"
                ) as f:
                    f.write(bytesio_object.getbuffer())
        else:
            return

    def _is_fname_bad(self, fname: str) -> bool:
        bad_fnames = self._load_bad_fnames()
        return fname in bad_fnames

    def _is_input_bad(self, *args) -> bool:
        efname = self._set_fnames(*args)
        return self._is_fname_bad(efname)

    def _load_bad_fnames(self) -> list[str]:
        with open(self._bad_input_path, "rb") as f:
            return pickle.load(f)

    def _store_bad_fname(self, fname: str) -> None:
        x = self._load_bad_fnames()
        x.append(fname)
        with open(self._bad_input_path, "wb") as f:
            pickle.dump(x, f)

    def _store_bad_input(self, *args) -> None:
        self._store_bad_fname(self._set_fnames(*args))

    def _load(self, *args) -> Union[np.ndarray, Dataset]:
        efname = self._set_fnames(*args)
        fpath = os.path.join(self._storage_dir, efname)
        if self.save_format == "hgt":
            with rasterio.open(fpath) as src:
                return src.read().squeeze()
        elif self.save_format == "nc":
            return Dataset(fpath, "r", format="NETCDF4")

    def _load_from_args(self, args: tuple[tuple], /, squeeze: bool = True):
        data = []
        for input_set in _strip_tuple(args):
            if self._is_input_bad(*input_set):
                print(f"Loading tile for {input_set} skipped, tile does not exist")
                data.append(None)
                continue
            data.append(self._load(*input_set))
        return data

    def _download_from_args(self, *args: Union[tuple[tuple], tuple]):
        for input_set in _strip_tuple(args):
            if self._before_request(*input_set):
                self._make_request(*input_set)
                self._after_request()