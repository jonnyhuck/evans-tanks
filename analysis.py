import os
os.environ['USE_PYGEOS'] = '0'

from affine import Affine
from geopandas import read_file
from numpy import sum as npsum
from rasterio.windows import Window
from rasterio import open as rio_open
from rasterio.features import rasterize

import warnings
warnings.filterwarnings("ignore")

# load tanks
tanks = read_file("./data/tanks_reprojected.shp")

# load MNDWI data
with rio_open("./data/NDWI_reproject.tiff") as ds:

    # loop through tanks
    results = []
    for id, tank in tanks.iterrows():

        # buffer
        buffer = tank.geometry.buffer(300)

        # get the top left & bottom right corner of ther AOI in image space
        bounds = buffer.bounds
        tl_img = ds.index(bounds[0], bounds[3])
        br_img = ds.index(bounds[2], bounds[1])
        w, h = br_img[1]-tl_img[1], br_img[0]-tl_img[0]

        # read using window for speed / RAM and create affine transform
        ndwi_band = ds.read(1, window=Window(tl_img[1], tl_img[0], w, h))
        affine = Affine(ds.res[0], 0, bounds[0], 0, -ds.res[1], bounds[3])

        # reclass
        ndwi_band[ndwi_band <= 0] = 0   # not wet
        ndwi_band[ndwi_band > 0] = 1    # wet

        # rasterize the buffer (we don't want to mask as it will involve operations on the lcm dataset directly)
        outer = rasterize([buffer], ndwi_band.shape, fill=0, transform=affine, default_value=1) 
        buffer_sum = npsum(outer * ndwi_band)

        # count of cells in tank
        inner = rasterize([tank.geometry], ndwi_band.shape, fill=0, transform=affine, default_value=1)
        tank_sum = npsum(inner)
        
        # work out the proportion (max 1) and store
        results.append(buffer_sum / tank_sum if tank_sum > buffer_sum else 1)

# load results into shapefile
tanks['evans_coeff'] = results

# export
tanks.to_file("./out/tank_wetness.shp")
print("done")