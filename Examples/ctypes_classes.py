import ctypes
import time
import numpy as np
from numpy.ctypeslib import ndpointer

# Create a structure for complex numbers
class Complex(ctypes.Structure):
   _fields_ = [("real", ctypes.c_double),
               ("imag", ctypes.c_double)
              ]

# load the C++ shared library as lib
lib_ccc = ctypes.cdll.LoadLibrary('../bin/ccc.so')

class CCC(object):
    def __init__(self, a, b, th, m2, m3, length):
        lib_ccc.ccc_new.argtypes = [ctypes.c_double,
                                    ctypes.c_double,
                                    ctypes.c_double,
                                    ctypes.c_double,
                                    ctypes.c_double,
                                    ctypes.c_int
                                   ]
        lib_ccc.ccc_new.restype = ctypes.c_void_p

        lib_ccc.get_cc.argtypes = [ctypes.c_void_p] 
        lib_ccc.get_cc.restypes =  ctypes.c_void_p

        lib_ccc.get_ca.argtypes = [ctypes.c_void_p]
        lib_ccc.get_ca.restypes =  ctypes.c_void_p

        lib_ccc.print_ccc.argtypes = [ctypes.c_void_p,ctypes.c_char_p]
        lib_ccc.print_ccc.restypes =  ctypes.c_void_p

        lib_ccc.copy_cc_ca.argtypes = [ctypes.c_void_p,
                                       ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                       ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS")]
        lib_ccc.copy_cc_ca.restypes = ctypes.c_void_p

        lib_ccc.copy_lenses.argtypes = [ctypes.c_void_p,
                                        ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS")
                                       ]
        lib_ccc.copy_lenses.restypes = ctypes.c_void_p

        lib_ccc.get_bounding_box.argtypes = [ctypes.c_void_p,
                                             ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                             ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                             ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                             ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                             ctypes.c_double
                                            ]
        lib_ccc.get_bounding_box.restypes = ctypes.c_void_p

        self.obj = lib_ccc.ccc_new(a, b, th, m2, m3, length)

    def get_cc(self):
        lib_ccc.get_cc(self.obj)

    def get_ca(self):
        lib_ccc.get_ca(self.obj)

    def print_ccc(self, file_name):
        lib_ccc.print_ccc(self.obj, file_name)

    def copy_cc_ca(self,ccp, cap):
        lib_ccc.copy_cc_ca(self.obj, ccp, cap)

    def copy_lenses(self, lens_pos):
        lib_ccc.copy_lenses(self.obj, lens_pos)

    def get_bounding_box(self, cc_min, cc_max, ca_min, ca_max, scale):
        lib_ccc.get_bounding_box(self.obj, cc_min, cc_max, ca_min, ca_max, scale)


# load the C++ shared library as lib
lib_img = ctypes.cdll.LoadLibrary('../bin/imgpoint.so')

# Represent functional ctypes version of C++ CCC class as Python Class
class Img(object):
    def __init__(self, a, b, th, m2, m3, pos_x, pos_y):
        lib_img.img_new.argtypes = [ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double
                               ]
        lib_img.img_new.restype = ctypes.c_void_p

        lib_img.get_roots.argtypes = [ctypes.c_void_p] 
        lib_img.get_roots.restypes = ctypes.c_void_p

        lib_img.get_images.argtypes = [ctypes.c_void_p]
        lib_img.get_images.restypes = ctypes.c_void_p

        lib_img.set_pos.argtypes = [ctypes.c_void_p,
                                ctypes.c_double,
                                ctypes.c_double
                               ]
        lib_img.set_pos.restype = ctypes.c_void_p

        lib_img.copy_images.argtypes = [ctypes.c_void_p,
                                    ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                    ndpointer(dtype=np.bool_,flags="C_CONTIGUOUS")]
        lib_img.copy_images.restypes = ctypes.c_void_p

        self.obj = lib_img.img_new(a, b, th, m2, m3, pos_x, pos_y)

    def get_roots(self):
        lib_img.get_roots(self.obj)

    def get_images(self):
        lib_img.get_images(self.obj)

    def set_pos(self, pos_x, pos_y):
        lib_img.set_pos(self.obj, pos_x, pos_y)

    def copy_images(self,roots, isimgs):
        lib_img.copy_images(self.obj, roots, isimgs)


# load the C++ shared library as lib
lib_lc = ctypes.cdll.LoadLibrary('../bin/lcirs.so')

# Represent functional ctypes version of C++ LightCurveIRS class as Python Class
class LC_irs(object):
    def __init__(self, a, b, th, m2, m3, source_size, lc_len, img_plane_size):
        lib_lc.lcirs_new.argtypes = [ctypes.c_double,
                                     ctypes.c_double,
                                     ctypes.c_double,
                                     ctypes.c_double,
                                     ctypes.c_double,
                                     ctypes.c_double,
                                     ctypes.c_uint,
                                     ctypes.c_long
                                    ]
        lib_lc.lcirs_new.restype = ctypes.c_void_p

        lib_lc.get_lc.argtypes = [ctypes.c_void_p,
                                  ctypes.c_double,
                                  ctypes.c_double,
                                  ctypes.c_double,
                                  ctypes.c_double
                                 ] 
        lib_lc.get_lc.restypes = ctypes.c_void_p

        lib_lc.get_lc_irs.argtypes = [ctypes.c_void_p,
                                      ctypes.c_double,
                                      ctypes.c_double,
                                      ctypes.c_double,
                                      ctypes.c_double
                                     ] 
        lib_lc.get_lc_irs.restypes = ctypes.c_void_p

        lib_lc.copy_lc.argtypes = [ctypes.c_void_p,
                                   ndpointer(dtype=np.double,flags="C_CONTIGUOUS")
                                  ]
        lib_lc.copy_lc.restypes = ctypes.c_void_p

        self.obj = lib_lc.lcirs_new(a, b, th, m2, m3, source_size, lc_len, img_plane_size)

    def get_lc(self, iniX, iniY, finX, finY):
        lib_lc.get_lc(self.obj, iniX, iniY, finX, finY)

    def get_lc_irs(self, iniX, iniY, finX, finY):
        lib_lc.get_lc_irs(self.obj, iniX, iniY, finX, finY)

    def copy_lc(self,lc_vec):
        lib_lc.copy_lc(self.obj, lc_vec)

