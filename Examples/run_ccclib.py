# This is a simple script to create an class based on the shared library.
# It creates a simple object with critical curve and caustic and then plots it as an png.

import ctypes
import time
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt

# load the C++ shared library as lib
lib = ctypes.cdll.LoadLibrary('../bin/ccc.so')

# Create a structure for complex numbers
class Complex(ctypes.Structure):
   _fields_ = [("real", ctypes.c_double),
               ("imag", ctypes.c_double)
              ]

# Represent functional ctypes version of C++ CCC class as Python Class
class CCC(object):
    def __init__(self, a, b, th, m2, m3, length):
        lib.ccc_new.argtypes = [ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_double,
                                ctypes.c_int
                               ]
        lib.ccc_new.restype = ctypes.c_void_p

        lib.get_cc.argtypes = [ctypes.c_void_p] 
        lib.get_cc.restypes = ctypes.c_void_p

        lib.get_ca.argtypes = [ctypes.c_void_p]
        lib.get_ca.restypes = ctypes.c_void_p

        lib.print_ccc.argtypes = [ctypes.c_void_p,ctypes.c_char_p]
        lib.print_ccc.restypes = ctypes.c_void_p

        self.obj = lib.ccc_new(a, b, th, m2, m3, length)

        lib.copy_cc_ca.argtypes = [ctypes.c_void_p,
                                   ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS"),
                                   ndpointer(dtype=np.cdouble,flags="C_CONTIGUOUS")]
        lib.copy_cc_ca.restypes = ctypes.c_void_p


    def get_cc(self):
        lib.get_cc(self.obj)

    def get_ca(self):
        lib.get_ca(self.obj)

    def print_ccc(self, file_name):
        lib.print_ccc(self.obj, file_name)

    def copy_cc_ca(self,ccp, cap):
        lib.copy_cc_ca(self.obj,ccp,cap)



# Define lens parameters.
a = 1.0
b = 1.0
theta = 1.05
m2 = 1/3
m3 = 1/3
length = 500

start_time = time.time()
# Initialise the CriticalCurveCaustic object.
ccc = CCC(a,b,theta,m2,m3,length)
ccc.get_ca()
ccc.print_ccc(ctypes.c_char_p(("./test.dat").encode('utf-8')))

# Copy the data-points from CCC object to numpy arrays
cc_array = np.zeros(3000, np.cdouble)
ca_array = np.zeros(3000, np.cdouble)
ccc.copy_cc_ca(cc_array, ca_array)

cc_real = []
cc_imag = []
ca_real = []
ca_imag = []

for cc in cc_array:
  cc_real.append(cc.real)
  cc_imag.append(cc.imag)

for ca in ca_array:
  ca_real.append(ca.real)
  ca_imag.append(ca.imag)

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.scatter(cc_real, cc_imag)
ax2.scatter(ca_real, ca_imag)
fig.savefig("CCC_test.png", dpi=100)

print("Time to initialise, calculate curves, copy & print & plot the data (s): ",time.time()-start_time)
