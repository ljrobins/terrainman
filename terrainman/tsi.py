
import os
import datetime
import requests
from bs4 import BeautifulSoup
import numpy as np
from scipy.interpolate import interp1d
from .util import DataProduct

class TsiDataHandler(DataProduct):
    """Handles Total Solar Irradiance (TSI) data

    :param DataProduct: Abstract class for handling data products
    :type DataProduct: class
    """
    def __init__(self) -> None:
        self._storage_dir = os.path.join(os.environ["TERRAIN_DATA"], "tsi")
        self.url_base = (
            "https://www.ncei.noaa.gov/data/total-solar-irradiance/access/daily/"
        )
        self.download_format = "nc"
        self.save_format = "nc"
        self._REF_DATE = datetime.datetime(
            year=1610,
            month=1,
            day=1,
            hour=0,
            minute=0,
            second=0,
            tzinfo=datetime.timezone.utc,
        )

        self._storage_file = os.path.join(os.environ["TERRAIN_DATA"], "tsi_daily.csv")
        if not os.path.exists(self._storage_file):
            page = requests.get(
                "https://www.ncei.noaa.gov/data/total-solar-irradiance/access/daily/"
            )
            soup = BeautifulSoup(page.content, "html.parser")
            # get all hrefs
            links = []
            for link in soup.find_all("a"):
                if "tsi_" in link.get("href"):
                    links.append(link.get("href"))

            with open(os.path.join(self._storage_file), "w") as f:
                f.write(",".join(links))

        self._after_init()

    def _set_fnames(self, date: datetime.datetime) -> str:
        csv_path = os.path.join(os.environ["TERRAIN_DATA"], f"tsi_daily.csv")
        with open(csv_path, "r") as f:
            csv_els = [x for x in f.read().split(",") if len(x) > 0]
            try:
                fnames = [c for c in csv_els if f"daily_s{date.year}" in c][-1]
            except IndexError:
                fnames = csv_els[-1]

        self.extracted_fname = fnames.split("/")[-1]
        self.url_fname = fnames
        return self.extracted_fname

    def download(self, date: datetime.datetime):
        """Downloads TSI data for a given date

        :param date: Date [utc]
        :type date: datetime.datetime
        """
        self._download_from_args((date,))

    def load(self, date: datetime.datetime):
        """Loads TSI data for a given date

        :param date: Date [utc]
        :type date: datetime.datetime
        :return: TSI data
        :rtype: _type_
        """
        return self._load_from_args((date,))[0]

    def _dates_to_total_seconds(self, dates: np.ndarray[datetime.datetime]):
        return np.array([(d - self._REF_DATE).total_seconds() for d in dates])

    def build_interp(self, years: np.ndarray[int]):
        """Builds an interpolation function for TSI data

        :param years: Unique years to build interpolation function for
        :type years: np.ndarray[int]
        """
        rel_day_vals = []
        tsi_vals = []
        for year in years:
            dset = self.load(datetime.datetime(year=year, month=1, day=1))

            tsi_vals = np.concatenate(
                (tsi_vals, dset.variables["TSI"][:].astype(np.float64).flatten())
            )

            tvals = dset.variables["time"][:].astype(np.float64)
            rel_day_vals = np.concatenate((rel_day_vals, tvals * 86400))

        self.interp = lambda dates: interp1d(
            rel_day_vals, tsi_vals, fill_value=1361.0, bounds_error=False
        )(self._dates_to_total_seconds(dates))

    def eval(self, dates: np.ndarray[datetime.datetime]) -> np.ndarray:
        """Evaluates TSI data for a given set of dates

        :param dates: Dates [utc]
        :type dates: np.ndarray[datetime.datetime]
        :return: Interpolated TSI data [W/m^2]
        :rtype: np.ndarray
        """
        uyears = np.unique([d.year for d in dates])
        [self.download(datetime.datetime(year=y, month=1, day=1)) for y in uyears]
        self.build_interp(uyears)
        return self.interp(dates)