# This is a simple script to create an class based on the shared library.
# It creates a simple object with critical curve and caustic and then plots it as an png.

import ctypes
import time
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from run_ccclib import CCC

# load the C++ shared library as lib
lib_lc = ctypes.cdll.LoadLibrary('../bin/lcirs.so')

# Create a structure for complex numbers
class Complex(ctypes.Structure):
   _fields_ = [("real", ctypes.c_double),
               ("imag", ctypes.c_double)
              ]

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
pos_x = 1.0
pos_y = 1.0

# number of steps
lc_steps = 100
source_size = 1e-3

lc_array     = np.zeros(lc_steps, np.double)
lc_array_irs = lc_array

lc_irs = LC_irs(a,b,theta, m2, m3, source_size, lc_steps, 300)
#lc_irs.get_lc(0.0,0.0,1.0,0.577)
#lc_irs.copy_lc(lc_array)

lc_irs.get_lc_irs(0.0,0.0,1.0,0.577)
lc_irs.copy_lc(lc_array)

print("Copied LC")
print(lc_array)

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.set_title("Source Trajectory")

ax1.axis(xmin=cc_min.real, xmax=cc_max.real, ymin=cc_min.imag, ymax=cc_max.imag)
ax1.scatter(ca_real, ca_imag, s = 0.1)
ax1.scatter(lenses_real, lenses_imag, s=200.0, marker = 'o')
ax1.plot([0.0,1.0],[0.0,0.577])

ax2.set_title("Light Curve")

ax2.plot(lc_array)

fig.savefig("LightCurvePoint.png", dpi=200)


#for i in range(0,100):
#  pos_x = 0.01*i
#  pos_y = 0.00577*i
#
#  print("position number ", str(i))  
#
#  img_generator.set_pos(pos_x, pos_y)
#  img_generator.get_images()
#  #img_generator = Img(a,b,theta,m2,m3,pos_x, pos_y)
#  #img_generator.get_images()
#
#  img_generator.copy_images(roots_array, isimg_array)
#
#  true_img_real = []
#  true_img_imag = []
#  fake_img_real = []
#  fake_img_imag = []
#  
#  for j in range(0,10):
#    if isimg_array[j]:
#      true_img_real.append(roots_array[j].real)
#      true_img_imag.append(roots_array[j].imag)
#    else:
#      fake_img_real.append(roots_array[j].real)
#      fake_img_imag.append(roots_array[j].imag)
#  
#  # Plotting
#  fig, (ax1, ax2) = plt.subplots(1, 2)
#  
#  title = "Number of true images: " + str(len(true_img_real))
#  
#  ax1.set_title(title)
#  
#  ax1.axis(xmin=cc_min.real, xmax=cc_max.real, ymin=cc_min.imag, ymax=cc_max.imag)
#  ax1.scatter(cc_real, cc_imag, s = 0.1)
#  ax1.scatter(lenses_real, lenses_imag, s=200.0, marker = 'o')
#  ax1.scatter(fake_img_real, fake_img_imag, s=100.0, color = 'red', marker = 'x')
#  ax1.scatter(true_img_real, true_img_imag, s=100.0, color = 'green', marker = 'x')
#  
#  pos_summary = ""
#  for k in range(0,len(true_img_real)):
#      pos_summary += "[" + f'{true_img_real[k]:.2f}' + ", " + f'{true_img_imag[k]:.2f}' + "]"
#      if (k % 3 == 2):
#          pos_summary += "\n"
#  
#  ax2.set_title(pos_summary, fontsize = 8)
#  
#  ax2.axis(xmin=ca_min.real, xmax=ca_max.real, ymin=ca_min.imag, ymax=ca_max.imag)
#  ax2.scatter(ca_real, ca_imag, s = 0.1)
#  ax2.scatter(lenses_real, lenses_imag, s = 200.0, marker = 'o')
#  ax2.scatter(pos_x, pos_y, s=100.0, color = 'green', marker = 'x')
#  
#  fig.savefig("Img_test"+str(i)+".png", dpi=200)
#
print("Time to initialise, calculate curves, copy & print & plot the data (s): ",time.time()-start_time)
