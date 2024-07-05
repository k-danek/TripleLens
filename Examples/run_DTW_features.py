import ctypes
from ctypes import *
import time
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from dtaidistance import dtw
from sklearn.cluster import KMeans
import math
import sys
from matplotlib.patches import Rectangle, Patch

from ctypes_classes import CCC
from ctypes_classes import LC_irs
from analytic_classes import LightCurveAnalyzer, Complex, get_trajectory, LightCurveDTW

sys.setrecursionlimit(500000) 


# Define lens parameters.
a = 1.01
b = 1.0001
theta = 1.047197551
m2 = 5.0e-3
m3 = 0.0
length = 500


# Define source parameters.
source_size = 1e-4
lc_steps = 512
points_per_radius = 30
ini_time = -1.0
fin_time = 1.0

q_values = [ 0.1, 0.15, 0.2, 0.3]  # Example values for q
#alpha_values = [np.pi / 8.0, np.pi / 5.0, np.pi / 2.0, 5.0 * np.pi / 8.0, 3.0 * np.pi / 4.0, 7.0 * np.pi / 8.0]  # Example values for alpha
alpha_values =  np.arange(0, 2 * np.pi, np.pi / 8.0)  # Example values for alpha

feature_names = ["caustic_entry", "caustic_exit", "cusp_approach_a", "cusp_approach_b", "cusp_transversal", "dip", "double_crossing"]

analyzerDTW = LightCurveDTW(a, b, theta, m2, m3, source_size, lc_steps, points_per_radius, feature_names)

# Color map for different short signals
#colors = plt.cm.get_cmap('tab10', len(feature_names))
#colors = plt.cm.tab10.colors

for q in q_values:
    for alpha in alpha_values:
        analyzerDTW.run_for_params(q, alpha, ini_time, fin_time)