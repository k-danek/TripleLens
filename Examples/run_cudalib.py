# This is a simple script to create an class based on the shared library.
# It creates a simple object with critical curve and caustic and then plots it as an png.

import ctypes
import time
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from ctypes_classes import CCC
from ctypes_classes import LC_cuda
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
pos_ini_x = 0.0
pos_ini_y = 0.0

# number of steps
lc_steps = 100
source_size = 1e-3
points_per_radius = 50

lc_point_array  = np.zeros(lc_steps, np.double)
lc_cuda_array   = np.zeros(lc_steps, np.double)
lc_irs_array    = np.zeros(lc_steps, np.double)
lc_irs_ac_array = np.zeros(lc_steps, np.double)

lc_cuda   = LC_cuda(a,b,theta, m2, m3, source_size, lc_steps, points_per_radius)
lc_irs    = LC_irs(a,b,theta, m2, m3, source_size, lc_steps, points_per_radius)
lc_irs_ac = LC_irs(a,b,theta, m2, m3, source_size, lc_steps, 4*points_per_radius)


temp_temp = time.time()
lc_cuda.get_lc_cuda(0.0,0.0,1.0,0.577)
lc_cuda.copy_lc(lc_cuda_array)
time_to_cuda = time.time() - temp_temp

temp_temp = time.time()
lc_irs.get_lc(0.0,0.0,1.0,0.577)
lc_irs.copy_lc(lc_point_array)
time_to_point = time.time() - temp_temp

temp_temp = time.time()
lc_irs.get_lc_irs(0.0,0.0,1.0,0.577)
lc_irs.copy_lc(lc_irs_array)
time_to_irs = time.time() - temp_temp


temp_temp = time.time()
lc_irs_ac.get_lc_irs(0.0,0.0,1.0,0.577)
lc_irs_ac.copy_lc(lc_irs_ac_array)
time_to_irs_ac = time.time() - temp_temp

print("Copied LC")

# Plotting
#fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 5)

ax1 = plt.subplot(221)
#ax1.set_title("Source Trajectory")
ax1.text(.5,.9,'Source Trajectory', horizontalalignment='center',transform=ax1.transAxes)
ax1.axis(xmin=cc_min.real, xmax=cc_max.real, ymin=cc_min.imag, ymax=cc_max.imag)
ax1.scatter(ca_real, ca_imag, s = 0.1)
ax1.scatter(lenses_real, lenses_imag, s=200.0, marker = 'o')
ax1.plot([0.0,1.0],[0.0,0.577])

ax2 = plt.subplot(422)
ax2.text(.7,.85,'Light Curve', horizontalalignment='center',transform=ax2.transAxes)
#ax2.set_title("Light Curve")

ax2.plot(lc_point_array, color='red', label='point')
ax2.plot(lc_irs_array, color='blue', label='irs')
ax2.plot(lc_cuda_array, color='green', label='cuda')

ax3 = plt.subplot(424)
ax3.text(.7,.85,'point vs irs', horizontalalignment='center',transform=ax3.transAxes)
ax3.plot(lc_point_array-lc_irs_array, color='purple', label='point')

ax4 = plt.subplot(426)
ax4.text(.7,.85,'point vs cuda', horizontalalignment='center',transform=ax4.transAxes)
ax4.plot(lc_point_array-lc_cuda_array, color='brown', label='point')

ax5 = plt.subplot(428)
ax5.text(.7,.85,'cuda vs irs', horizontalalignment='center',transform=ax5.transAxes)
ax5.plot(lc_cuda_array-lc_irs_array, color='cyan', label='point')

ax6 = plt.subplot(425)
ax6.text(.7,.85,'precise irs vs cuda', horizontalalignment='center',transform=ax6.transAxes)
ax6.plot(lc_irs_ac_array-lc_cuda_array, color='black', label='point')

ax7 = plt.subplot(427)
ax7.text(.7,.85,'precise irs vs irs', horizontalalignment='center',transform=ax7.transAxes)
ax7.plot(lc_irs_ac_array-lc_irs_array, color='black', label='point')



plt.savefig("LightCurvePointCuda.png", dpi=400)

print("Time to initialise and calculate the images (s): ",time.time()-start_time)
print("Time to compute cuda irs (s): ",time_to_cuda)
print("Time to compute standard irs (s): ",time_to_irs)
print("Time to compute higher accuracy irs (s): ",time_to_irs_ac)
print("Time to compute a point amplification (s): ",time_to_point)
