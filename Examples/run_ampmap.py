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
import math
from decimal import Decimal, ROUND_HALF_UP

# Create a structure for complex numbers
class Complex(ctypes.Structure):
   _fields_ = [("real", ctypes.c_double),
               ("imag", ctypes.c_double)
              ]


def grid_to_coordinate_(min, max, grid_size, grid_point):
  return (min+(max-min)*float(grid_point)/float(grid_size)) 

def coordinate_to_grid(min, max, grid_size, coordinate_point):
  return int(float(grid_size-1)*(coordinate_point-min)/(max-min)) 

# Define lens parameters, wide
a = 5.0
b = 5.0
theta = 1.047197551
m2 = 1.0e-1
m3 = 1.0e-1
length = 500


start_time = time.time()
# Initialise the CriticalCurveCaustic object.
ccc = CCC(a,b,theta,m2,m3,length)
ccc.get_ca()

# copy the bounding box
cc_min = np.zeros(1, np.cdouble)
cc_max = np.zeros(1, np.cdouble)
ca_min = np.zeros(1, np.cdouble)
ca_max = np.zeros(1, np.cdouble)
ccc.get_bounding_box(cc_min, cc_max, ca_min, ca_max, 1.2)
print("Bounding box cc: ", cc_min, cc_max)
print("Bounding box ca: ", ca_min, ca_max)


# imgs

pos_ini_x = ca_min.real
pos_fin_x = ca_max.real
pos_ini_y = ca_min.imag
pos_fin_y = ca_max.imag


# number of steps
lc_steps = 500
source_size = 1e-3
points_per_radius = 30

lc_point_array = np.zeros(lc_steps, np.double)
lc_irs_array   = np.zeros(lc_steps, np.double)

lc_irs = LC_irs(a, b, theta, m2, m3, source_size, lc_steps, points_per_radius)

data_array = np.empty((0, lc_steps), dtype=np.float64)

for i in range(0, lc_steps):
  row_y = pos_ini_y+float(i/(lc_steps-1.0))*(pos_fin_y-pos_ini_y)
  lc_irs.get_lc(pos_ini_x,row_y,pos_fin_x,row_y)
  lc_irs.copy_lc(lc_point_array)
  data_array = np.vstack([data_array, lc_point_array])


print("Bounding box cc: ", cc_min, cc_max)
print("Bounding box ca: ", ca_min, ca_max)

print("bottom left corner",coordinate_to_grid(pos_ini_x, pos_fin_x, lc_steps, pos_ini_x), coordinate_to_grid(pos_ini_y, pos_fin_y, lc_steps, pos_ini_y))
print("origin",coordinate_to_grid(pos_ini_x, pos_fin_x, lc_steps, 0.0), coordinate_to_grid(pos_ini_y, pos_fin_y, lc_steps, 0.0))
print("top right corner",coordinate_to_grid(pos_ini_x, pos_fin_x, lc_steps, pos_fin_x), coordinate_to_grid(pos_ini_y, pos_fin_y, lc_steps, pos_fin_y))


plt.rcParams.update({'font.size': 8}) 

fig, axs = plt.subplots(1, 2)

v_min = 1.0
v_max = 10.0
# Plotting the heatmap
img = axs[1].imshow(data_array, cmap='plasma', vmin=v_min, vmax=v_max, aspect='auto')


lc_ini = -0.5
lc_fin = 5.5
lc_ini_x = coordinate_to_grid(pos_ini_x, pos_fin_x, lc_steps, lc_ini)
lc_fin_x = coordinate_to_grid(pos_ini_x, pos_fin_x, lc_steps, lc_fin)
lc_ini_y = coordinate_to_grid(pos_ini_y, pos_fin_y, lc_steps, 0.0)
lc_fin_y = coordinate_to_grid(pos_ini_y, pos_fin_y, lc_steps, 0.0)
plt.plot([lc_ini_x, lc_fin_x], [lc_ini_y, lc_fin_y], color="green", linewidth=2)


# Plotting source points
lens_1_x = coordinate_to_grid(pos_ini_x, pos_fin_x, lc_steps, 0.0)
lens_1_y = coordinate_to_grid(pos_ini_y, pos_fin_y, lc_steps, 0.0)
lens_2_x = coordinate_to_grid(pos_ini_x, pos_fin_x, lc_steps, a)
lens_2_y = coordinate_to_grid(pos_ini_y, pos_fin_y, lc_steps, 0.0)
lens_3_x = coordinate_to_grid(pos_ini_x, pos_fin_x, lc_steps, b*math.cos(theta))
lens_3_y = coordinate_to_grid(pos_ini_y, pos_fin_y, lc_steps, b*math.sin(theta))
axs[1].scatter(lens_1_x, lens_1_y, marker='o', facecolors='none', edgecolors='white', s=100)
axs[1].scatter(lens_2_x, lens_2_y, marker='o', facecolors='none', edgecolors='white', s=100)
axs[1].scatter(lens_3_x, lens_3_y, marker='o', facecolors='none', edgecolors='white', s=100)


axs[1].set_aspect('equal', adjustable='box')

# Ascenting y indices are ascending on y axis now
axs[1].invert_yaxis()

# Add colorbar for reference
cbar = plt.colorbar(img, fraction=0.045, orientation='vertical', pad=0.05)
cbar.mappable.set_clim(v_min, v_max)

# Set your custom x and y axis ranges
label_step = 1.0

pos_ini_x_rounded = float(np.round(pos_ini_x, 1))
pos_ini_y_rounded = float(np.round(pos_ini_y, 1))

x_label_range = np.arange(pos_ini_x_rounded, pos_fin_x, label_step)
y_label_range = np.arange(pos_ini_y_rounded, pos_fin_y, label_step)

x_grid_range = []
y_grid_range = []

for i in range(0,len(x_label_range)):
  x_grid_range.append(coordinate_to_grid(pos_ini_x, pos_fin_x, lc_steps, x_label_range[i]))

for i in range(0,len(y_label_range)):
  y_grid_range.append(coordinate_to_grid(pos_ini_y, pos_fin_y, lc_steps, y_label_range[i]))

print("lenghth of grid and labels", len(x_grid_range), len(x_label_range))

axs[1].set_xticks(x_grid_range)
axs[1].set_xticklabels([f'{float(val):.1f}' for val in x_label_range])
axs[1].set_yticks(y_grid_range)
axs[1].set_yticklabels([f'{float(val):.1f}' for val in y_label_range])

# Get the lc 
lc_irs.get_lc_irs(lc_ini,-0.00,lc_fin,-0.00)
lc_irs.copy_lc(lc_point_array)

lc_range_array = np.arange(lc_ini, lc_fin, (lc_fin-lc_ini)/float(lc_steps))
axs[0].plot(lc_range_array, lc_point_array, color="green", linewidth=2)
axs[0].set_box_aspect(0.905)

fig.savefig('amplification_map.svg', format='svg', bbox_inches='tight')
fig.savefig('amplification_map.pdf', format='pdf', bbox_inches='tight')
fig.savefig('amplification_map.eps', format='eps', bbox_inches='tight')


# Show the plot
plt.show() 

print("Copied LC")