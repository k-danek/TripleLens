# This is a simple script to create an class based on the shared library.
# It creates a simple object with critical curve and caustic and then plots it as an png.

import ctypes
import time
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from ctypes_classes import CCC

# Create a structure for complex numbers
class Complex(ctypes.Structure):
   _fields_ = [("real", ctypes.c_double),
               ("imag", ctypes.c_double)
              ]

# Define lens parameters.
a = 1.5
b = 1.5
theta = 1.05
m2 = 1/3
m3 = 1/3
length = 500

start_time = time.time()
# Initialise the CriticalCurveCaustic object.
ccc = CCC(a,b,theta,m2,m3,length)
ccc.get_ca()

# Print out Critical Curve & Caustics
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

# copy lens positions
lens_pos = np.zeros(3, np.cdouble)
ccc.copy_lenses(lens_pos)
lenses_real = []
lenses_imag = []
for lens in lens_pos:
  lenses_real.append(lens.real)
  lenses_imag.append(lens.imag)

# copy the bounding box
cc_min = np.zeros(1, np.cdouble)
cc_max = np.zeros(1, np.cdouble)
ca_min = np.zeros(1, np.cdouble)
ca_max = np.zeros(1, np.cdouble)
ccc.get_bounding_box(cc_min, cc_max, ca_min, ca_max, 1.1)
print("Bounding box cc: ", cc_min, cc_max)
print("Bounding box ca: ", ca_min, ca_max)

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.axis(xmin=cc_min.real, xmax=cc_max.real, ymin=cc_min.imag, ymax=cc_max.imag)
ax1.scatter(cc_real, cc_imag)
ax1.scatter(lenses_real, lenses_imag, s= 200.0, marker = 'o')

ax2.axis(xmin=ca_min.real, xmax=ca_max.real, ymin=ca_min.imag, ymax=ca_max.imag)
ax2.scatter(ca_real, ca_imag)
ax2.scatter(lenses_real, lenses_imag, s = 200.0, marker = 'o')

fig.savefig("CCC_test.png", dpi=100)

print("Time to initialise, calculate curves, copy & print & plot the data (s): ",time.time()-start_time)
