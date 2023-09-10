# This is a simple script to create an class based on the shared library.
# It creates a simple object with critical curve and caustic and then plots it as an png.

import ctypes
import time
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from ctypes_classes import CCC
from ctypes_classes import LC_irs
import csv
import scipy.cluster.hierarchy as hcluster
from sklearn.cluster import MeanShift

class LineSegment:
  def __init__(self, nY, nL, nR):
    self.nY = nY
    self.nL = nL
    self.nR = nR


lineSegments = []


with open('Amoeba_150', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in spamreader:
        lineSegments.append(LineSegment(int(row[0]),int(row[1]),int(row[2])))

print("read csv")

num_of_segments = len(lineSegments)
boundaries = np.ndarray(shape=(2*len(lineSegments),2), dtype=int, order='F')

for i in range(0,num_of_segments):
  #np.concatenate((boundaries,np.array([[int(row[0])],[int(row[1])]], dtype=int)), axis=1)
  #np.concatenate((boundaries,np.array([[int(row[0])],[int(row[2])]], dtype=int)), axis=1)
  boundaries[2*i][0]   = lineSegments[i].nY
  boundaries[2*i][1]   = lineSegments[i].nL
  boundaries[2*i+1][0] = lineSegments[i].nY
  boundaries[2*i+1][1] = lineSegments[i].nR


mean_shift = MeanShift()
mean_shift.fit(boundaries)
cluster_labels = mean_shift.labels_
num_of_clusters = len(np.unique(cluster_labels))

print("Cluster labels " + str(cluster_labels))
print("Num of clusters " + str(num_of_clusters))

# clustering
#thresh = 10
#clusters = hcluster.fclusterdata(boundaries, thresh, criterion="distance")

#num_of_clusters = max(clusters)
#print("num of clusters "+str(num_of_clusters))
#print("len of clusters "+str(len(clusters)))
#print("vals of clusters "+str(clusters[0])+","+str(clusters[1])+","+str(clusters[2]))

#print("lineSegments:"+str(lineSegments[0].nY)+","+str(lineSegments[0].nL))
#print("boundaries:"+str(boundaries[0][0])+","+str(boundaries[0][1]))
#print("boundaries:"+str(boundaries[1][0])+","+str(boundaries[1][1]))
#print("boundaries:"+str(boundaries[2][0])+","+str(boundaries[2][1]))
#print("boundaries lenght:"+str(boundaries.size))
#print("boundaries shape:"+str(boundaries.shape))

#clustered_boundaries = [[] * (num_of_clusters)]
#clustered_boundaries = scipy.zeros((n,1), 'double')
clustered_boundaries = []

# Adding spurious elements because of inferior knowledge of Python.
for i in range(0,num_of_clusters):
  clustered_boundaries.append([LineSegment(0,0,0)])

for i in range(0,len(lineSegments)):
  clustered_boundaries[cluster_labels[2*i]].append(LineSegment(boundaries[2*i][0], boundaries[2*i][1], boundaries[2*i+1][1]))




## Adding spurious elements because of inferior knowledge of Python.
#for i in range(0,num_of_clusters):
#  clustered_boundaries.append([LineSegment(0,0,0)])
#
##print("boundaries shape:"+str(clustered_boundaries))
#
#for i in range(0,len(lineSegments)):
#  clustered_boundaries[clusters[i]-1].append(LineSegment(boundaries[2*i][0], boundaries[2*i][1], boundaries[2*i+1][1]))
#
#clusters_to_delete = []

#for i in range(0,num_of_clusters):
#  if len(clustered_boundaries[i]) > 1:
#    clustered_boundaries[i].pop(0)
#    print("clustered boundary "+str(i)+" length: "+str(len(clustered_boundaries[i])))
#  else:
#    clusters_to_delete.append(i)
#
#for i in clusters_to_delete:
#    del clustered_boundaries[i]
index = 0
while index < len(clustered_boundaries):
  if len(clustered_boundaries[index]) > 1:
    clustered_boundaries[index].pop(0)
    print("clustered boundary "+str(index)+" length: "+str(len(clustered_boundaries[index])))
    index = index+1
  else:
    print("clustered boundary "+str(index)+" length: "+str(len(clustered_boundaries[index]))+" deleted")
    del clustered_boundaries[index]  


print("Clustered boundaries:"+str(len(clustered_boundaries)))
#print("Clustered boundaries:"+str(clustered_boundaries))

# Plotting
#fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(2, 4)

ax1 =plt.subplot(221)
ax1.set_title("Image Overview")
for segment in lineSegments:
  ax1.plot([segment.nL,segment.nR],[segment.nY,segment.nY])

ax2 = plt.subplot(422)
for segment in clustered_boundaries[0]:
  ax2.plot([segment.nL,segment.nR],[segment.nY,segment.nY])

if len(clustered_boundaries) > 1:
  ax3 = plt.subplot(424)
  for segment in clustered_boundaries[1]:
    ax3.plot([segment.nL,segment.nR],[segment.nY,segment.nY])

if len(clustered_boundaries) > 2:
  ax4 = plt.subplot(426)
  for segment in clustered_boundaries[2]:
    ax4.plot([segment.nL,segment.nR],[segment.nY,segment.nY])

if len(clustered_boundaries) > 3:
  ax5 = plt.subplot(428)
  for segment in clustered_boundaries[3]:
    ax5.plot([segment.nL,segment.nR],[segment.nY,segment.nY])

if len(clustered_boundaries) > 4:
  ax6 = plt.subplot(425)
  for segment in clustered_boundaries[4]:
    ax6.plot([segment.nL,segment.nR],[segment.nY,segment.nY])

if len(clustered_boundaries) > 5:
  ax7 = plt.subplot(427)
  for segment in clustered_boundaries[5]:
    ax7.plot([segment.nL,segment.nR],[segment.nY,segment.nY])


#ax2.set_title("Light Curve")
#
#ax2.plot(lc_point_array, color='cyan', label='point')
#ax2.plot(lc_irs_array, color='red', label='irs')
#
#ax3.set_title("Ratio Curve")
#
#ax3.plot(lc_irs_array/lc_point_array, color='cyan', label='point')
##ax3.set_ylim([1.0, 1.5])
plt.savefig("Amoeba.png", dpi=300)

