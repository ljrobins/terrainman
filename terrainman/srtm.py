import os
from typing import Union, Tuple
from .util import DataProduct
import numpy as np
import scipy


class TerrainDataHandler(DataProduct):
    """Handles Shuttle Radar Terrain Mission (SRTM) data

    :param DataProduct: Abstract class for handling data products
    :type DataProduct: class
    """

    def __init__(self) -> None:
        self._storage_dir = os.path.join(os.environ["TERRAIN_DATA"], "terrain")
        self.url_base = "https://e4ftl01.cr.usgs.gov/MEASURES/SRTMGL1.003/2000.02.11/"
        self.download_format = "zip"
        self.save_format = "hgt"
        self.HTTP_ERR_MSG = (
            "HTTP error occured, likely because tile is ocean, aborting..."
        )
        self.DNE_MSG = "Tile does not exist, likely because it's ocean, skipping"
        self.EXISTS_MSG = "Tile already downloaded, skipping"

        self._after_init()

    def _set_fnames(self, lat: int, lon: int) -> str:
        tile_lon_str = f'{"W" if lon < 0 else "E"}{abs(lon):03}'
        tile_lat_str = f'{"N" if lat > 0 else "S"}{abs(lat):02}'
        extracted_fname = f"{tile_lat_str}{tile_lon_str}.hgt"
        url_fname = f"{tile_lat_str}{tile_lon_str}.SRTMGL1.hgt.zip"
        self.extracted_fname = extracted_fname
        self.url_fname = url_fname
        return self.extracted_fname

    def download_tile(self, lat: int, lon: int) -> None:
        """Downloads [1 deg x 1 deg] a SRTM 30 meter resolution terrain tile for a given (integer) latitude and longitude pair, requires `os.environ['EARTHDATA_USERNAME']` and `os.environ['EARTHDATA_PASSWORD']` to be defined

        :param lat: Latitude [deg]
        :type lat: int
        :param lon: Longitude [deg]
        :type lon: int
        """
        self._download_from_args((lat, lon))

    def load_tile(self, lat: int, lon: int) -> np.ndarray:
        """Loads [1 deg x 1 deg] SRTM 30 meter resolution terrain tile for a given (integer) latitude and longitude pair

        :param lat: Geodetic latitude [deg]
        :type lat: int
        :param lon: Longitude [deg]
        :type lon: int
        :return: Elevation grid [m] above the geoid (MSL)
        :rtype: np.ndarray
        """
        elev_grid = self._load_from_args((lat, lon))[0]
        if elev_grid is not None:
            lat_space = (lat + 1) - np.linspace(0, 1, elev_grid.shape[0])
            lon_space = lon + np.linspace(0, 1, elev_grid.shape[0])
            return SrtmTile(lat_space, lon_space, elev_grid)
        return None

    def load_tiles_containing(
        self, lat: np.ndarray, lon: np.ndarray
    ) -> list[np.ndarray]:
        """Loads [1 deg x 1 deg] SRTM 30 meter resolution terrain tiles to cover all latitude and longitudes input

        :param lat: Geodetic latitudes [deg]
        :type lat: np.ndarray
        :param lon: Longitudes [deg]
        :type lon: np.ndarray
        :return: Elevation grids [m] above the geoid (MSL)
        :rtype: list[np.ndarray]
        """
        lat, lon = np.array([lat]).flatten(), np.array([lon]).flatten()
        self.download_tiles_containing(lat, lon)
        unique_lat_lons = self._unique_tile_lat_lon_pairs(lat, lon)
        return (
            [self.load_tile(*args) for args in unique_lat_lons]
            if len(unique_lat_lons) > 1
            else self.load_tile(*unique_lat_lons[0])
        )

    def download_tiles_containing(self, lat: np.ndarray, lon: np.ndarray) -> None:
        """Downloads [1 deg x 1 deg] SRTM 30 meter resolution terrain tiles to cover all latitude and longitudes input, requires `os.environ['EARTHDATA_USERNAME']` and `os.environ['EARTHDATA_PASSWORD']` to be defined

        :param lat: Latitude [deg]
        :type lat: np.ndarray
        :param lon: Longitude [deg]
        :type lon: np.ndarray
        """
        unique_lat_lons = self._unique_tile_lat_lon_pairs(lat, lon)
        self._download_from_args(unique_lat_lons)

    def _unique_tile_lat_lon_pairs(
        self, lat: np.ndarray, lon: np.ndarray
    ) -> list[Tuple[int, int]]:
        lat_lon_dec = np.vstack((lat.flatten(), lon.flatten())).T
        lat_lon_int = np.floor(lat_lon_dec)
        unique_lat_lons = np.unique(lat_lon_int, axis=0)
        unique_lat_lons = np.array(unique_lat_lons, np.int32)
        return tuple(map(tuple, unique_lat_lons))


class SrtmTile:
    """Represents a 1 x 1 degree terrain tile for Shuttle Radar Terrain Mission (SRTM) data"""

    def __init__(
        self, lat_space: np.ndarray, lon_space: np.ndarray, elev_grid: np.ndarray
    ) -> None:
        self.lat_space, self.lon_space = lat_space, lon_space
        self.lat_grid, self.lon_grid = np.meshgrid(lat_space, lon_space, indexing="ij")
        self.elev_grid = elev_grid

    def interpolate(
        self, lat_deg: Union[float, np.ndarray], lon_deg: Union[float, np.ndarray]
    ) -> Union[float, np.ndarray]:
        """Interpolates within the latitude/longitude

        :param lat_deg: Geodetic latitude [deg]
        :type lat_deg: Union[float, np.ndarray]
        :param lon_deg: Longitude [deg]
        :type lon_deg: Union[float, np.ndarray]
        :raises NotImplemented: Inputs must be either floats or np.ndarrays
        :return: Elevation [m] above the geoid (MSL)
        :rtype: Union[float, np.ndarray]
        """
        kwargs = {
            "points": (self.lat_space, self.lon_space),
            "values": self.elev_grid,
            "fill_value": np.nan,
            "bounds_error": False,
        }
        if isinstance(lat_deg, np.ndarray):
            return scipy.interpolate.interpn(
                **kwargs,
                xi=np.hstack(
                    (lat_deg.reshape(lat_deg.size, 1), lon_deg.reshape(lat_deg.size, 1))
                ),
            ).reshape(lat_deg.shape)
        elif isinstance(lat_deg, float):
            return scipy.interpolate.interpn(
                **kwargs, xi=np.hstack((lat_deg, lon_deg))
            )[0]
        else:
            raise NotImplementedError("Inputs must be either floats or np.ndarrays")
