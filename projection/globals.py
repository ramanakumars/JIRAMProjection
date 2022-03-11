import numpy as np
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spiceypy as spice
import netCDF4 as nc
import json, glob, re, os
import multiprocessing, time
from scipy.interpolate import interp2d, griddata
import ctypes

## load the C library to get the projection mask
print(os.path.dirname(__file__))
project_c = np.ctypeslib.load_library('project.so', os.path.dirname(__file__))

image_mask_c = project_c.get_image_mask
process_c = project_c.process

array_1d_int    = np.ctypeslib.ndpointer(dtype=np.int, ndim=1, flags='C_CONTIGUOUS')
array_1d_double = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='C_CONTIGUOUS')
array_2d_double = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='C_CONTIGUOUS')

image_mask_c.argtypes = [array_1d_double, array_1d_double, ctypes.c_int, ctypes.c_int,\
                         array_2d_double, array_2d_double, array_1d_double, ctypes.c_int]
image_mask_c.restype  = array_1d_int

process_c.argtypes = [ctypes.c_double, ctypes.c_int, array_1d_double, \
                         array_2d_double, array_2d_double, array_2d_double, array_2d_double, array_2d_double]


## and the spice furnish function for the library
furnish_c    = project_c.furnish
furnish_c.argtypes = [ctypes.c_char_p]

FRAME_HEIGHT = 128
IMG_HEIGHT   = 128
FRAME_WIDTH  = 432
PIXEL_SCALE  = 0.00023778 ## rad/pixel

FOCAL_LENGTH = 159.82033 ## m
PIXEL_SIZE   = 0.038 ## mm 
CENTER_X     = 216.5
CENTER_Y     = 64.5
F1           = FOCAL_LENGTH/PIXEL_SIZE

