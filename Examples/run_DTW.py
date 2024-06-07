from dtaidistance import dtw
from sklearn.cluster import KMeans
import numpy as np
import math

import matplotlib.pyplot as plt
from ctypes_classes import CCC
from ctypes_classes import LC_irs

def getTrajectory(q, alpha, iniTime, finTime, steps):
  trajectory = []

  for i in range (0, steps):
    time = (finTime-iniTime)*i/(steps-1)
    point = Complex(cos(alpha)*time-q*sin(alpha), sin(alpha)*time+q*cos(alpha))
    trajectory.append(point)

  return trajectory

# Define lens parameters.
a = 1.0
b = 1.0
theta = 1.047197551
m2 = 5.0e-2
m3 = 1.0e-3

# number of steps
lc_steps = 512
source_size = 1e-5
points_per_radius = 10

ini_time = -1.0
fin_time =  1.0


qs = np.arange(0.1, 0.31, 0.2)
alphas = np.arange(0, 6 * np.pi, np.pi / 4.0+0.01)
data = []

lc_irs = LC_irs(a,b,theta, m2, m3, source_size, lc_steps, points_per_radius)

for q in qs:
    for alpha in alphas:
        pos_ini_x = math.cos(alpha)*ini_time-q*math.sin(alpha)
        pos_ini_y = math.sin(alpha)*ini_time+q*math.cos(alpha)
        pos_fin_x = math.cos(alpha)*fin_time-q*math.sin(alpha)
        pos_fin_y = math.sin(alpha)*fin_time+q*math.cos(alpha)
        lc_irs.get_lc(pos_ini_x,pos_ini_y,pos_fin_x,pos_fin_y)
        lc_point_array = np.zeros(lc_steps, np.double)
        lc_irs.copy_lc(lc_point_array)
        data.append(lc_point_array)


print("data 0:"+str(data[0]))
print("data 10:"+str(data[10]))

def calculate_derivative(time_series):
    return np.gradient(time_series)

# Compute derivatives for all time series
data_derivative = [np.gradient(series) for series in data]

# Assuming `data` is your 1D time series dataset, shape (n_samples, 128)
distances = np.zeros((len(data), len(data)))
for i in range(len(data)):
    for j in range(i, len(data)):
        distances[i, j] = distances[j, i] = dtw.distance(data_derivative[i], data_derivative[j])

# Print distances in ascending order
all_distances = []



# Gather all distances with their indices
for i in range(len(data)):
    for j in range(i + 1, len(data)):  # Only consider upper triangle to avoid duplicates
        all_distances.append((distances[i, j], i, j))

# Sort all distances by value
all_distances.sort()

# Print sorted distances
print("\nSorted distances (value, i, j):")
for distance, i, j in all_distances:
    print(f"{distance:.4f} ({i},{j})")

# Set the backend to 'Agg'
plt.switch_backend('Agg')

# Function to plot 9 pairs
def plot_9_pairs(pairs, plot_num):
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    for k, (dist, i, j) in enumerate(pairs):
        ax = axes[k // 3, k % 3]
        path = dtw.warping_path(data[i], data[j])
        aligned_data_i = [data[i][index_i] for index_i, index_j in path]
        aligned_data_j = [data[j][index_j] for index_i, index_j in path]
        ax.plot(aligned_data_i, label=f"Series {i}")
        ax.plot(aligned_data_j, label=f"Series {j} (aligned)")
        ax.set_title(f"Pair {i}, {j} (Distance: {dist:.4f})")
        ax.legend()
    plt.tight_layout()
    plt.savefig(f"top_9_pairs_{plot_num}.png")
    plt.close()

# Plot all pairs in chunks of 9
num_pairs = len(all_distances)
#for plot_num in range((num_pairs + 8) // 27):  # +8 to ensure we cover all pairs
for plot_num in range((num_pairs + 8) // 9):  # +8 to ensure we cover all pairs
    start_index = plot_num * 9
    end_index = min(start_index + 9, num_pairs)
    plot_9_pairs(all_distances[start_index:end_index], plot_num)