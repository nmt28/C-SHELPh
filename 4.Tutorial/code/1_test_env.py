try:
    import numpy as np
    import h5py as h5
    from matplotlib.widgets import LassoSelector
    from matplotlib.path import Path
    import matplotlib.pyplot as plt
    import pyproj
    from pyproj import Proj
    from pyproj import Transformer
    import pandas as pd
    import argparse
    import os
    import subprocess
    import time
    import utm
    import math
    import fiona
    import geopandas
    import netCDF4
    from datetime import datetime

    print("Environment Test Successful")
except Exception as e:
 raise e
