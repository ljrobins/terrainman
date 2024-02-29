import os
import dotenv
import inspect

from .tile_io import *

if os.path.exists('.env.shared'):
    dotenv.load_dotenv('.env.shared')


os.environ['TERRAIN_SRC'] = os.path.dirname(
    os.path.abspath(inspect.getsourcefile(lambda: 0))
)
os.environ['TERRAIN_DATA'] = os.path.join(os.environ['TERRAIN_SRC'], 'data')
_DATA_PRODUCT_TYPES = ["aerosol", "cloud", "terrain", 'tsi', 'viirs_color']

if not os.path.exists(os.environ['TERRAIN_DATA']):
    os.mkdir(os.environ['TERRAIN_DATA'])
for product in _DATA_PRODUCT_TYPES:
    product_data_path = os.path.join(os.environ['TERRAIN_DATA'], product)
    if not os.path.exists(product_data_path):
        os.mkdir(product_data_path)