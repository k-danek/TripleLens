# This is a simple script to create an class based on the shared library.
# It creates a simple object with critical curve and caustic and then plots it as an png.

import ctypes
from ctypes import *
import time
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from ctypes_classes import CCC
from ctypes_classes import LC_irs

# Create a structure for complex numbers
class Complex(ctypes.Structure):
   _fields_ = [("real", ctypes.c_double),
               ("imag", ctypes.c_double)
              ]

# Define lens parameters.
a = 1.0
b = 1.0
theta = 1.047197551
m2 = 1/3
m3 = 1/3
length = 500

start_time = time.time()
# Initialise the CriticalCurveCaustic object.
ccc = CCC(a,b,theta,m2,m3,length)
ccc.get_ca()

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
ccc.get_bounding_box(cc_min, cc_max, ca_min, ca_max, 1.5)
print("Bounding box cc: ", cc_min, cc_max)
print("Bounding box ca: ", ca_min, ca_max)

# imgs
#pos_ini_x = 0.0
#pos_ini_y = 0.0
#pos_fin_x = 1.0
#pos_fin_y = 0.577

# imgs
pos_ini_x = 0.0
pos_ini_y = -0.7
pos_fin_x = 1.0
pos_fin_y = 0.7

# source pos
#pos_ini_x = 0.22
#pos_ini_y = 0.12694
#pos_fin_x = 0.24
#pos_fin_y = 0.13848

# number of steps
lc_steps = 200
source_size = 8e-3
points_per_radius = 100

lc_point_array = np.zeros(lc_steps, np.double)
lc_irs_array   = np.zeros(lc_steps, np.double)

lc_irs = LC_irs(a,b,theta, m2, m3, source_size, lc_steps, points_per_radius)

amoeba_filename_buffer = create_string_buffer(b"Amoeba_")
param_filename_buffer = create_string_buffer(b"Pars.dat")

lc_irs.set_amoeba_printout(amoeba_filename_buffer,param_filename_buffer)

lc_irs.get_lc_irs(pos_ini_x,pos_ini_y,pos_fin_x,pos_fin_y)
lc_irs.copy_lc(lc_irs_array)

lc_irs.get_lc(pos_ini_x,pos_ini_y,pos_fin_x,pos_fin_y)
lc_irs.copy_lc(lc_point_array)

print("Copied LC")

# Plotting
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

ax1.set_title("Source Trajectory")

ax1.axis(xmin=cc_min.real, xmax=cc_max.real, ymin=cc_min.imag, ymax=cc_max.imag)
ax1.scatter(ca_real, ca_imag, s = 0.1)
ax1.scatter(lenses_real, lenses_imag, s=200.0, marker = 'o')
ax1.plot([pos_ini_x,pos_fin_x],[pos_ini_y,pos_fin_y])

ax2.set_title("Light Curve")

ax2.plot(lc_point_array, color='cyan', label='point')
ax2.plot(lc_irs_array, color='red', label='irs')

ax3.set_title("Ratio Curve")

ax3.plot(lc_irs_array/lc_point_array, color='cyan', label='point')
#ax3.set_ylim([1.0, 1.5])
fig.savefig("LightCurvePoint.png", dpi=300)

print("Time to initialise and calculate the images (s): ",time.time()-start_time)
