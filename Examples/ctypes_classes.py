import ctypes
import time
import numpy as np
from numpy.ctypeslib import ndpointer
import copy

# Create a structure for complex numbers
class Complex(ctypes.Structure):
   _fields_ = [("real", ctypes.c_double),
               ("imag", ctypes.c_double)
              ]

class CCC(object):
    def __init__(self, a, b, th, m2, m3, length):
        # load the C++ shared library as lib
        self.lib_ccc = ctypes.cdll.LoadLibrary('../bin/libccc.so')

        self.lib_ccc.ccc_new.argtypes = [ctypes.c_double,
                                         ctypes.c_double,
                                         ctypes.c_double,
                                         ctypes.c_double,
                                         ctypes.c_double,
                                         ctypes.c_int
                                        ]
        self.lib_ccc.ccc_new.restype = ctypes.c_void_p

        self.lib_ccc.get_cc.argtypes = [ctypes.c_void_p] 
        self.lib_ccc.get_cc.restypes =  ctypes.c_void_p

        self.lib_ccc.get_ca.argtypes = [ctypes.c_void_p]
        self.lib_ccc.get_ca.restypes =  ctypes.c_void_p

        self.lib_ccc.print_ccc.argtypes = [ctypes.c_void_p,ctypes.c_char_p]
        self.lib_ccc.print_ccc.restypes =  ctypes.c_void_p

        self.lib_ccc.copy_cc_ca.argtypes = [ctypes.c_void_p,
                                            ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                            ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS")
                                           ]
        self.lib_ccc.copy_cc_ca.restypes = ctypes.c_void_p

        self.lib_ccc.copy_lenses.argtypes = [ctypes.c_void_p,
                                             ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS")
                                            ]
        self.lib_ccc.copy_lenses.restypes = ctypes.c_void_p

        self.lib_ccc.get_bounding_box.argtypes = [ctypes.c_void_p,
                                                  ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                                  ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                                  ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                                  ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                                  ctypes.c_double
                                                 ]
        self.lib_ccc.get_bounding_box.restypes = ctypes.c_void_p

        self.obj = self.lib_ccc.ccc_new(a, b, th, m2, m3, length)

    def get_cc(self):
        self.lib_ccc.get_cc(self.obj)

    def get_ca(self):
        self.lib_ccc.get_ca(self.obj)

    def print_ccc(self, file_name):
        self.lib_ccc.print_ccc(self.obj, file_name)

    def copy_cc_ca(self,ccp, cap):
        self.lib_ccc.copy_cc_ca(self.obj, ccp, cap)

    def copy_lenses(self, lens_pos):
        self.lib_ccc.copy_lenses(self.obj, lens_pos)

    def get_bounding_box(self, cc_min, cc_max, ca_min, ca_max, scale):
        self.lib_ccc.get_bounding_box(self.obj, cc_min, cc_max, ca_min, ca_max, scale)



# Represent functional ctypes version of C++ CCC class as Python Class
class Img(object):
    def __init__(self, a, b, th, m2, m3, pos_x, pos_y):
        # load the C++ shared library as lib
        self.lib_img = ctypes.cdll.LoadLibrary('../bin/libimg.so')

        self.lib_img.img_new.argtypes = [ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double
                               ]
        self.lib_img.img_new.restype = ctypes.c_void_p

        self.lib_img.get_roots.argtypes = [ctypes.c_void_p] 
        self.lib_img.get_roots.restypes = ctypes.c_void_p

        self.lib_img.get_images.argtypes = [ctypes.c_void_p]
        self.lib_img.get_images.restypes = ctypes.c_void_p

        self.lib_img.set_pos.argtypes = [ctypes.c_void_p,
                                         ctypes.c_double,
                                         ctypes.c_double
                                        ]
        self.lib_img.set_pos.restype = ctypes.c_void_p

        self.lib_img.copy_images.argtypes = [ctypes.c_void_p,
                                             ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                             ndpointer(dtype=np.bool_,flags="C_CONTIGUOUS")
                                            ]
        self.lib_img.copy_images.restypes = ctypes.c_void_p

        self.obj = self.lib_img.img_new(a, b, th, m2, m3, pos_x, pos_y)

    def get_roots(self):
        self.lib_img.get_roots(self.obj)

    def get_images(self):
        self.lib_img.get_images(self.obj)

    def set_pos(self, pos_x, pos_y):
        self.lib_img.set_pos(self.obj, pos_x, pos_y)

    def copy_images(self,roots, isimgs):
        self.lib_img.copy_images(self.obj, roots, isimgs)


# Represent functional ctypes version of C++ LightCurveIRS class as Python Class
class LC_irs(object):
    def __init__(self, a, b, th, m2, m3, source_size, lc_len, img_plane_size):
        # load the C++ shared library as lib
        self.lib_lc = ctypes.cdll.LoadLibrary('../bin/liblcirs.so')

        self.lib_lc.lcirs_new.argtypes = [ctypes.c_double,
                                          ctypes.c_double,
                                          ctypes.c_double,
                                          ctypes.c_double,
                                          ctypes.c_double,
                                          ctypes.c_double,
                                          ctypes.c_uint,
                                          ctypes.c_long
                                         ]
        self.lib_lc.lcirs_new.restype = ctypes.c_void_p

        self.lib_lc.get_lc.argtypes = [ctypes.c_void_p,
                                       ctypes.c_double,
                                       ctypes.c_double,
                                       ctypes.c_double,
                                       ctypes.c_double
                                      ] 
        self.lib_lc.get_lc.restypes = ctypes.c_void_p

        self.lib_lc.get_lc_irs.argtypes = [ctypes.c_void_p,
                                           ctypes.c_double,
                                           ctypes.c_double,
                                           ctypes.c_double,
                                           ctypes.c_double
                                          ] 
        self.lib_lc.get_lc_irs.restypes = ctypes.c_void_p

        self.lib_lc.copy_lc.argtypes = [ctypes.c_void_p,
                                        ndpointer(dtype=np.double,flags="C_CONTIGUOUS")
                                       ]
        self.lib_lc.copy_lc.restypes = ctypes.c_void_p

        self.lib_lc.set_limb_darkening.argtypes = [ctypes.c_void_p,
                                                   ctypes.c_char_p,
                                                   ctypes.c_double
                                                  ]
        self.lib_lc.set_limb_darkening.restypes = ctypes.c_void_p

        self.lib_lc.set_amoeba_printout.argtypes = [ctypes.c_void_p,
                                                    ctypes.c_char_p,
                                                    ctypes.c_char_p
                                                  ]
        self.lib_lc.set_amoeba_printout.restypes = ctypes.c_void_p


        self.obj = self.lib_lc.lcirs_new(a, b, th, m2, m3, source_size, lc_len, img_plane_size)

    def get_lc(self, iniX, iniY, finX, finY):
        self.lib_lc.get_lc(self.obj, iniX, iniY, finX, finY)

    def get_lc_irs(self, iniX, iniY, finX, finY):
        self.lib_lc.get_lc_irs(self.obj, iniX, iniY, finX, finY)

    def copy_lc(self,lc_vec):
        self.lib_lc.copy_lc(self.obj, lc_vec)

    def set_limb_darkening(self, model_type, model_parameter):
        self.lib_lc.set_limb_darkening(self.obj, model_type, model_parameter)

    def set_amoeba_printout(self, amoeba_filename, par_filename):
        self.lib_lc.set_amoeba_printout(self.obj, amoeba_filename, par_filename)


# Represent functional ctypes version of C++ LightCurveCUDA class as Python Class
class LC_cuda(object):
    def __init__(self, a, b, th, m2, m3, source_size, lc_len, img_plane_size):
        print("Python Init: called")
        
        # load the C++ shared library as lib
        self.lib_cuda = ctypes.cdll.LoadLibrary('../bin/liblccuda.so')

        print("Python Init: linked binary")

        self.lib_cuda.lccuda_new.argtypes = [ctypes.c_double,
                                             ctypes.c_double,
                                             ctypes.c_double,
                                             ctypes.c_double,
                                             ctypes.c_double,
                                             ctypes.c_double,
                                             ctypes.c_uint,
                                             ctypes.c_long
                                            ]
        self.lib_cuda.lccuda_new.restype = ctypes.c_void_p

        self.lib_cuda.get_lc.argtypes = [ctypes.c_void_p,
                                         ctypes.c_double,
                                         ctypes.c_double,
                                         ctypes.c_double,
                                         ctypes.c_double
                                        ] 
        self.lib_cuda.get_lc.restypes = ctypes.c_void_p


        self.lib_cuda.get_lc_cuda.argtypes = [ctypes.c_void_p,
                                              ctypes.c_double,
                                              ctypes.c_double,
                                              ctypes.c_double,
                                              ctypes.c_double
                                             ] 
        self.lib_cuda.get_lc_cuda.restypes = ctypes.c_void_p

        self.lib_cuda.copy_lc.argtypes = [ctypes.c_void_p,
                                          ndpointer(dtype=np.double,flags="C_CONTIGUOUS")
                                         ]
        self.lib_cuda.copy_lc.restypes = ctypes.c_void_p

        self.lib_cuda.set_limb_darkening_cuda.argtypes = [ctypes.c_void_p,
                                                          ctypes.c_char_p,
                                                          ctypes.c_double
                                                         ]
        self.lib_cuda.set_limb_darkening_cuda.restypes = ctypes.c_void_p

        print("Python Init: defined all the functions")

        self.obj = self.lib_cuda.lccuda_new(copy.copy(a),
                                            copy.copy(b),
                                            copy.copy(th),
                                            copy.copy(m2),
                                            copy.copy(m3),
                                            copy.copy(source_size),
                                            copy.copy(lc_len),
                                            copy.copy(img_plane_size)
                                            )

        print("Python Init: called lccuda_new")


    def __del__(self):
        print("LC_Cuda destroyed")

    def get_lc(self, iniX, iniY, finX, finY):
        self.lib_cuda.get_lc(self.obj, iniX, iniY, finX, finY)

    def get_lc_cuda(self, iniX, iniY, finX, finY):
        self.lib_cuda.get_lc_cuda(self.obj, iniX, iniY, finX, finY)

    def copy_lc(self,lc_vec):
        self.lib_cuda.copy_lc(self.obj, lc_vec)

    def set_limb_darkening(self, model_type, model_parameter):
        self.lib_cuda.set_limb_darkening_cuda(self.obj, model_type, model_parameter)
