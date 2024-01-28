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
from ctypes_classes import LC_cuda

# Create a structure for complex numbers
class Complex(ctypes.Structure):
   _fields_ = [("real", ctypes.c_double),
               ("imag", ctypes.c_double)
              ]

def getCCC(a,b,theta,m2,m3,length):
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

  return ccc, cc_real, cc_imag, ca_real, ca_imag


# Define lens parameters.
a = 1.0
b = 1.0
theta = 1.047197551
m2 = 1/3
m3 = 1/30
length = 500


# Initialise the CriticalCurveCaustic object.

ccc_triple, cc_real, cc_imag, ca_real, ca_imag = getCCC(a,b,theta,m2,m3,length)
_, cc_real_double, cc_imag_double, ca_real_double, ca_imag_double = getCCC(a,b,theta,m2,0.0,length)

# copy the bounding box
cc_min = np.zeros(1, np.cdouble)
cc_max = np.zeros(1, np.cdouble)
ca_min = np.zeros(1, np.cdouble)
ca_max = np.zeros(1, np.cdouble)
ccc_triple.get_bounding_box(cc_min, cc_max, ca_min, ca_max, 1.5)
print("Bounding box cc: ", cc_min, cc_max)
print("Bounding box ca: ", ca_min, ca_max)

# copy lens positions
lens_pos = np.zeros(3, np.cdouble)
ccc_triple.copy_lenses(lens_pos)
lenses_real = []
lenses_imag = []
for lens in lens_pos:
  lenses_real.append(lens.real)
  lenses_imag.append(lens.imag)

# imgs
pos_ini_x = 0.0
pos_ini_y = -0.7
pos_fin_x = 1.0
pos_fin_y = 0.7

# number of steps
lc_steps = 200
source_size = 8e-3
points_per_radius = 100

lc_point_array_double = np.zeros(lc_steps, np.double)
lc_irs_array_double   = np.zeros(lc_steps, np.double)

lc_point_array_triple = np.zeros(lc_steps, np.double)
lc_irs_array_triple   = np.zeros(lc_steps, np.double)

# Change to LC_cuda if full CUDA testing is needed
lc_irs_triple = LC_irs(a,b,theta, m2, m3, source_size, lc_steps, points_per_radius)
lc_irs_double = LC_irs(a,b,theta, m2, 0.0, source_size, lc_steps, points_per_radius)

start_time = time.time()
lc_irs_triple.get_lc_irs(pos_ini_x,pos_ini_y,pos_fin_x,pos_fin_y)
lc_irs_triple.copy_lc(lc_irs_array_triple)
time_triple_irs = time.time()-start_time

start_time = time.time()
lc_irs_triple.get_lc(pos_ini_x,pos_ini_y,pos_fin_x,pos_fin_y)
lc_irs_triple.copy_lc(lc_point_array_triple)
time_triple_point = time.time()-start_time

start_time = time.time()
lc_irs_double.get_lc_irs(pos_ini_x,pos_ini_y,pos_fin_x,pos_fin_y)
lc_irs_double.copy_lc(lc_irs_array_double)
time_double_irs = time.time()-start_time

start_time = time.time()
lc_irs_double.get_lc(pos_ini_x,pos_ini_y,pos_fin_x,pos_fin_y)
lc_irs_double.copy_lc(lc_point_array_double)
time_double_point = time.time()-start_time

print("Copied LC")

# Plotting
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

ax1.set_title("Source Trajectory")

ax1.axis(xmin=cc_min.real, xmax=cc_max.real, ymin=cc_min.imag, ymax=cc_max.imag)
ax1.scatter(lenses_real, lenses_imag, s=200.0, marker = 'o')
ax1.scatter(ca_real_double, ca_imag_double, s = 0.2)
ax1.scatter(ca_real, ca_imag, s = 0.1)
ax1.plot([pos_ini_x,pos_fin_x],[pos_ini_y,pos_fin_y])


ax2.set_title("Light Curve")

ax2.plot(lc_point_array_triple, color='cyan', label='point')
ax2.plot(lc_irs_array_triple, color='red', label='irs')

ax2.plot(lc_point_array_double, color='green', label='point_double')
ax2.plot(lc_irs_array_double, color='blue', label='irs_double')

ax3.set_title("Ratio Curve")

ax3.plot((lc_irs_array_triple-lc_point_array_triple)/lc_point_array_triple, color='cyan', label='point')
ax3.plot((lc_irs_array_double-lc_point_array_double)/lc_point_array_double, color='green', label='point')
fig.savefig("LightCurveBinary.png", dpi=300)

print("Time to get light curve for triple lens with point source (s): ", time_triple_point)
print("Time to get light curve for triple lens with extended source (s): ", time_triple_irs)
print("Time to get light curve for double lens with point source (s): ", time_double_point)
print("Time to get light curve for double lens with extended source (s): ", time_double_irs)
