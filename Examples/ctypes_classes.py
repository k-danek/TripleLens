import ctypes
import time
import numpy as np
from numpy.ctypeslib import ndpointer

# load the C++ shared library as lib
lib = ctypes.cdll.LoadLibrary('../bin/imgpoint.so')

# Create a structure for complex numbers
class Complex(ctypes.Structure):
   _fields_ = [("real", ctypes.c_double),
               ("imag", ctypes.c_double)
              ]

# Represent functional ctypes version of C++ CCC class as Python Class
class Img(object):
    def __init__(self, a, b, th, m2, m3, pos_x, pos_y):
        lib.img_new.argtypes = [ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double
                               ]
        lib.img_new.restype = ctypes.c_void_p

        lib.get_roots.argtypes = [ctypes.c_void_p] 
        lib.get_roots.restypes = ctypes.c_void_p

        lib.get_images.argtypes = [ctypes.c_void_p]
        lib.get_images.restypes = ctypes.c_void_p

        lib.set_pos.argtypes = [ctypes.c_void_p,
                                ctypes.c_double,
                                ctypes.c_double
                               ]
        lib.set_pos.restype = ctypes.c_void_p

        lib.copy_images.argtypes = [ctypes.c_void_p,
                                    ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                    ndpointer(dtype=np.bool_,flags="C_CONTIGUOUS")]
        lib.copy_images.restypes = ctypes.c_void_p

        self.obj = lib.img_new(a, b, th, m2, m3, pos_x, pos_y)

    def get_roots(self):
        lib.get_roots(self.obj)

    def get_images(self):
        lib.get_images(self.obj)

    def set_pos(self, pos_x, pos_y):
        lib.set_pos(self.obj, pos_x, pos_y)

    def copy_images(self,roots, isimgs):
        lib.copy_images(self.obj, roots, isimgs)

