# This is a simple script to create an class based on the shared library.
# It creates a simple object with critical curve and caustic and then plots it as an png.

import ctypes
from ctypes import *
import time
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
import csv
import scipy.cluster.hierarchy as hcluster
from sklearn.cluster import MeanShift
from operator import attrgetter
import json
import matplotlib.transforms

from ctypes_classes import CCC

pars_file = "Par"
amoeba_file = "Amoeba_"

class LineSegment:
  def __init__(self, nY, nL, nR):
    self.nY = nY
    self.nL = nL
    self.nR = nR

class ParManager:
  def __init__(self, filename):
    par_dict = self.read_par_file(filename)
    #par_tuple = self.read_par_file(filename)
    #self.bottom_left_x = par_tuple[0]
    #self.bottom_left_y = par_tuple[1]
    #self.grid_step = par_tuple[2]
    self.bottom_left_x = float(par_dict['bottomLeftCornerImgX'])
    self.bottom_left_y = float(par_dict['bottomLeftCornerImgY'])
    self.grid_step = float(par_dict['imgPlaneSizeDouble'])/float(par_dict['imgPlaneSize'])
    self.z1_x = float(par_dict['z1x'])
    self.z1_y = float(par_dict['z1y'])
    self.z2_x = float(par_dict['z2x'])
    self.z2_y = float(par_dict['z2y'])
    self.z3_x = float(par_dict['z3x'])
    self.z3_y = float(par_dict['z3y'])

  def read_par_file(self,filename):
    file = open(filename, "r")
    parameters =  json.load(file)
    return parameters

  def ny_to_y(self, ny):
    return self.bottom_left_y + ny*self.grid_step
  def nx_to_x(self, nx):
    return self.bottom_left_x + nx*self.grid_step


def Get_bounding_box_segments(segments, enlarging_coeff):
  min_x = min(segments,key=attrgetter('nL'))
  max_x = max(segments,key=attrgetter('nR'))
  min_y = min(segments,key=attrgetter('nY'))
  max_y = max(segments,key=attrgetter('nY'))

  half_range_x = (max_x.nR-min_x.nL)/2
  half_range_y = (max_y.nY-min_y.nY)/2
  half_range_uni = enlarging_coeff * max(half_range_x, half_range_y)
  center_x = (max_x.nR+min_x.nL)/2
  center_y = (max_y.nY+min_y.nY)/2
  return [center_x-half_range_uni, center_x+half_range_uni, center_y-half_range_uni, center_y+half_range_uni,]

# Define lens parameters.
a = 1.22
b = 1.22
theta = 1.047197551
m2 = 1/3
m3 = 1/3
length = 1500

# imgs
pos_ini_x = 0.61
pos_ini_y = -0.1
pos_fin_x = 0.61
pos_fin_y = 0.85

# number of steps
lc_steps = 200
source_size = 1e-2
points_per_radius = 50



source_pos = []
for i in range(0, lc_steps):
  pos_x = pos_ini_x+(pos_fin_x-pos_ini_x)*float(i)/float(lc_steps-1)
  pos_y = pos_ini_y+(pos_fin_y-pos_ini_y)*float(i)/float(lc_steps-1)
  source_pos.append([pos_x, pos_y])

# Initialise the CriticalCurveCaustic object.
ccc = CCC(a,b,theta,m2,m3,length)
ccc.get_ca()

# Print out Critical Curve & Caustics
ccc.print_ccc(ctypes.c_char_p(("./test.dat").encode('utf-8')))

# Copy the data-points from CCC object to numpy arrays
cc_array       = np.zeros(6*length, np.cdouble)
ca_array       = np.zeros(6*length, np.cdouble)
ccc.copy_cc_ca(cc_array, ca_array)

print("copied ccc:"+str(len(cc_array)))
print("copied ccc:"+str(cc_array[0].real))
print("copied ccc:"+str(cc_array[0].imag))

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

par_manager = ParManager("Pars_Art.dat")

source_index = 79

lineSegments = []
index = 0
filename = "Amoeba_"+"Article_"+str(source_index)
with open(filename, newline='') as csvfile:
    print("file opened")
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    print("streamrider")
    for row in spamreader:
        index = index +1
        lineSegments.append(LineSegment(par_manager.ny_to_y(int(row[0])),
                                        par_manager.nx_to_x(int(row[1])),
                                        par_manager.nx_to_x(int(row[2]))))
        print("row read:" + str(index) )

print("read csv")

num_of_segments = len(lineSegments)

print("num of segments:" + str(num_of_segments))
boundaries = np.ndarray(shape=(2*len(lineSegments),2), dtype=float, order='F')

print("Defined boundaries")

for i in range(0,num_of_segments):
  boundaries[2*i][0]   = lineSegments[i].nY
  boundaries[2*i][1]   = lineSegments[i].nL
  boundaries[2*i+1][0] = lineSegments[i].nY
  boundaries[2*i+1][1] = lineSegments[i].nR

print("Copied boundaries")

clustered_boundaries = []

#mean_shift = MeanShift(bandwidth=0.1)
#mean_shift.fit(boundaries)
#cluster_labels = mean_shift.labels_
#num_of_clusters = len(np.unique(cluster_labels))
#
#clustered_boundaries = []
#
## Adding spurious elements because of inferior knowledge of Python.
#for i in range(0,num_of_clusters):
#  clustered_boundaries.append([LineSegment(0,0,0)])
#
#for i in range(0,len(lineSegments)):
#  clustered_boundaries[cluster_labels[2*i]].append(LineSegment(boundaries[2*i][0], boundaries[2*i][1], boundaries[2*i+1][1]))

# DO NOT DELETE - an image blob production.
#for cluster_index in range(0, len(clustered_boundaries)):
#  with open('clustered_boundaries_'+str(source_index)+'_'+str(cluster_index), 'w') as f:
#    for segment in clustered_boundaries[cluster_index]:
#      f.write(str(segment.nY)+','+str(segment.nR)+','+str(segment.nL)+'\n')
#    f.close()

#index = 0
#while index < len(clustered_boundaries):
#  if len(clustered_boundaries[index]) > 1:
#    clustered_boundaries[index].pop(0)
#    print("clustered boundary "+str(index)+" length: "+str(len(clustered_boundaries[index])))
#    index = index+1
#  else:
#    print("clustered boundary "+str(index)+" length: "+str(len(clustered_boundaries[index]))+" deleted")

# clustered_boundaries_79_6
for cluster_index in range(0,7):
  filename = "clustered_boundaries_79_"+str(cluster_index)
  with open(filename, newline='') as csvfile:
    # Adding spurious elements because of inferior knowledge of Python.
    clustered_boundaries.append([LineSegment(0,0,0)])
    clustered_boundaries[cluster_index].pop(0)
    print("file opened")
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    print("streamrider")
    for row in spamreader:
      clustered_boundaries[cluster_index].append(LineSegment(
                                                    float(row[0]),
                                                    float(row[1]),
                                                    float(row[2])))
    clustered_boundaries[cluster_index].pop(0)

print("Clustered boundaries:"+str(len(clustered_boundaries)))

# Plotting
#subplot_list = [(2,1),(3,1),(0,2),(1,2),(2,2),(3,2),(0,3),(1,3),(2,3),(3,3)]
subplot_list = [(1,2),(0,3),(1,3)]
axis_list = []
#color_list = ["blue", "green", "red", "cyan", "magenta", "yellow", '#0080ff', '#aa0511', '#0f0f0f', '#123456'] 
color_list = ["blue","blue", "blue","blue","blue","blue","blue","blue"] 
plt.rcParams.update({'font.size': 8, 'font.family': 'font.arial', 'svg.fonttype': 'none'})

#https://stackoverflow.com/questions/14600948/matplotlib-plot-outputs-text-as-paths-and-cannot-be-converted-to-latex-by-inks
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['figure.figsize'] = [12, 5]


fig = plt.figure()
dx = 2/72.; dy = 2/72.
scale_trans = fig.dpi_scale_trans
x_tic_offset = matplotlib.transforms.ScaledTranslation(0, dy, scale_trans)
y_tic_offset = matplotlib.transforms.ScaledTranslation(dx, 0, scale_trans)
ax1 = plt.subplot(121)
ax1.set_aspect('equal')
ax1.scatter(ca_real, ca_imag, color='red', s=0.01)
ax1.scatter(cc_real, cc_imag, color='blue', s=0.01)
ax1.scatter(source_pos[source_index][0], source_pos[source_index][1], marker='x', color='black', s=50)
ax1.scatter(par_manager.z1_x, par_manager.z1_y, marker='+', color='black', s=50)
ax1.scatter(par_manager.z2_x, par_manager.z2_y, marker='+', color='black', s=50)
ax1.scatter(par_manager.z3_x, par_manager.z3_y, marker='+', color='black', s=50)
#https://matplotlib.org/3.1.1/tutorials/intermediate/gridspec.html
for cb_index in range(0, len(clustered_boundaries)):
  for segment in clustered_boundaries[cb_index]:
    ax1.plot([segment.nL,segment.nR],[segment.nY,segment.nY], color="green")
for tick in ax1.xaxis.get_majorticklabels():
  tick.set_transform(tick.get_transform() + x_tic_offset)
for tick in ax1.yaxis.get_majorticklabels():
  tick.set_transform(tick.get_transform() + y_tic_offset)

# plotting detail of the source
ax2 = plt.subplot2grid((2,4), (0,2))
dx = 2/72.; dy = 2/72. 
x_tic_offset = matplotlib.transforms.ScaledTranslation(0, dy, scale_trans)
y_tic_offset = matplotlib.transforms.ScaledTranslation(dx, 0, scale_trans)
ax2.scatter(ca_real, ca_imag, color='red', s=0.02)
#ax2.plot(ca_real, ca_imag, color='#000033')
source_circle = plt.Circle((source_pos[source_index][0], source_pos[source_index][1]), source_size, color='black')
ax2.add_patch(source_circle)
source_detail_halfsize = source_size * 3
ax2.set_aspect('equal')
ax2.set_xlim([source_pos[source_index][0]-source_detail_halfsize, source_pos[source_index][0]+source_detail_halfsize])
ax2.set_ylim([source_pos[source_index][1]-source_detail_halfsize, source_pos[source_index][1]+source_detail_halfsize]) 
for tick in ax2.xaxis.get_majorticklabels():
  tick.set_transform(tick.get_transform() + x_tic_offset)
for tick in ax2.yaxis.get_majorticklabels():
  tick.set_transform(tick.get_transform() + y_tic_offset)

##axis_list.append(plt.subplot(subplot_list[i]))
for i in range(0,3):
  axis_list.append(plt.subplot2grid((2,4), subplot_list[i]))

  for tick in axis_list[i].xaxis.get_majorticklabels():
    tick.set_verticalalignment("top")
  bounding_box = Get_bounding_box_segments(clustered_boundaries[i], 1.15)
  for segment in clustered_boundaries[i]:
    axis_list[i].plot([segment.nL,segment.nR],[segment.nY,segment.nY], color='green')
  axis_list[i].scatter(cc_real, cc_imag, color='blue', s=0.02)
  axis_list[i].set_aspect('equal')
  axis_list[i].set_xlim([bounding_box[0], bounding_box[1]])
  axis_list[i].set_ylim([bounding_box[2], bounding_box[3]])
  for tick in axis_list[i].xaxis.get_majorticklabels():
    tick.set_transform(tick.get_transform() + x_tic_offset)
  for tick in axis_list[i].yaxis.get_majorticklabels():
    tick.set_transform(tick.get_transform() + y_tic_offset)



plt.savefig(filename+".pdf", dpi=600)
plt.clf()
plt.cla()
plt.close()

