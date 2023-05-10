import os
import dotenv

from .tile_io import *

if os.path.exists('.env.shared'):
    dotenv.load_dotenv('.env.shared')