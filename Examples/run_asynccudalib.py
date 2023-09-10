# This is a simple script to create an class based on the shared library.
# It creates a simple object with critical curve and caustic and then plots it as an png.

import ctypes
import time
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from ctypes_classes import CCC
from ctypes_classes import LC_cuda
#from ctypes_classes import LC_irs
from multiprocessing import Process

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


num_of_cores = 2

pos_ini_x = 0.0
pos_ini_y = 0.0
pos_fin_x = 1.0
pos_fin_y = 0.577

src_ini_x = [0.0] * num_of_cores
src_ini_y = [0.0] * num_of_cores
src_fin_x = [0.0] * num_of_cores
src_fin_y = [0.0] * num_of_cores

lc_cuda_instances = [LC_cuda] * num_of_cores

# number of steps
lc_steps = 100
lc_sub_steps = int(lc_steps / num_of_cores)
source_size = 1e-3
points_per_radius = 30

def run_cuda_job(lc_cuda_instance, ini_x,ini_y,fin_x,fin_y, cuda_out_array):
    lc_cuda_instance.get_lc_cuda(ini_x,ini_y,fin_x,fin_y)
    lc_cuda_instance.copy_lc(cuda_out_array)

#lc_cuda_array   = np.zeros(lc_steps, np.double)
lc_cuda_array   = []
lc_cuda_arrays  = np.zeros(shape=[lc_sub_steps]*num_of_cores, dtype=np.double)
process_array = [Process] * num_of_cores


temp_temp = time.time()

for i in range(0,num_of_cores):
    src_ini_x[i] = pos_ini_x + (pos_fin_x - pos_ini_x)*i/num_of_cores
    src_fin_x[i] = pos_ini_x + (pos_fin_x - pos_ini_x)*(i+1)/num_of_cores
    src_ini_y[i] = pos_ini_x + (pos_fin_y - pos_ini_y)*i/num_of_cores
    src_fin_y[i] = pos_ini_y + (pos_fin_y - pos_ini_y)*(i+1)/num_of_cores
    lc_cuda_instances[i] = LC_cuda(a,b,theta, m2, m3, source_size, lc_sub_steps, points_per_radius)
    process_array[i] = Process(target=run_cuda_job, args=(lc_cuda_instances[i],
                                                          src_ini_x[i],
                                                          src_fin_x[i],
                                                          src_ini_y[i],
                                                          src_fin_y[i],
                                                          lc_cuda_arrays[i]))
    process_array[i].start()


for i in range(0,num_of_cores):
    process_array[i].join()
    lc_cuda_array.append(lc_cuda_arrays[i])

time_to_cuda = time.time() - temp_temp

print("Time to compute cuda irs (s): ",time_to_cuda)


## Plotting
##fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 5)
#fig, (ax1, ax2) = plt.subplots(1, 2)
#
#ax1.set_title("Source Trajectory")
#
#ax1.axis(xmin=cc_min.real, xmax=cc_max.real, ymin=cc_min.imag, ymax=cc_max.imag)
#ax1.scatter(ca_real, ca_imag, s = 0.1)
#ax1.scatter(lenses_real, lenses_imag, s=200.0, marker = 'o')
#ax1.plot([0.0,1.0],[0.0,0.577])
#
#ax2.set_title("Light Curve")
#
#ax2.plot(lc_cuda_array)
#
#fig.savefig("LightCurvePoint.png", dpi=200)
#
#plt.savefig("LightCurvePointCuda.png", dpi=400)

