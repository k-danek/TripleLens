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
import csv
import json


from ctypes_classes import CCC
from ctypes_classes import LC_irs
from analytic_classes import LightCurveAnalyzer, Complex, get_trajectory, LightCurveDTW

sys.setrecursionlimit(1500000) 


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


# Define parameter ranges.
#a_values = [0.9, 0.95, 1.0, 1.1]  # Example values for a
#q_values = [0.05, 0.1, 0.15, 0.2, 0.4]  # Example values for q
#alpha_values = np.arange(0, 2 * np.pi, np.pi / 4.0)  # Example values for alpha

a_values = [0.90]  # Example values for a
q_values = [0.05, 0.1, 0.15, 0.2, 0.4]  # Example values for q
#q_values = [0.05]  # Example values for q
alpha_values = np.arange(0, 2 * np.pi, np.pi / 16.0) # Example values for alpha
#alpha_values = [np.pi / 2.0] # Example values for alpha


#feature_names = ["caustic_entry", "caustic_exit", "cusp_approach_a", "cusp_approach_b", "cusp_transversal", "dip", "double_crossing"]

# Output files
csv_file = "feature_matrix.csv"
json_file = "parameters_all.json"
feature_file_name = "features_summary"

# JSON structure to hold parameter values and line ranges
json_data = []

start_line = 0
end_line = 0


# broken ones : q=0.05, alpha = 5.49779, 3.927, 3.141529 
# super bronken: q=0.1, alpha= 0.785, 2.356

with open("feature_matrix", 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    for a in a_values:
        analyzerDTW = LightCurveDTW(a, b, theta, m2, m3, source_size, lc_steps, points_per_radius, feature_file_name)

        for q in q_values:
            for alpha in alpha_values:
                [feature_vector, light_curve, residual] = analyzerDTW.run_for_params(q, alpha, ini_time, fin_time, end_line % 1 == 0)
                writer.writerow(feature_vector)
                end_line += 1
                json_data.append({
                    "a": a,
                    "b": b,
                    "theta": theta,
                    "m2": m2,
                    "m3": m3,
                    "q": q,
                    "alpha": alpha,
                    "feature_vector": feature_vector,
                    "light_curve": residual.tolist()
                })

        #json_data.append({
        #            "a": a,
        #            "b": b,
        #            "theta": theta,
        #            "m2": m2,
        #            "m3": m3,
        #            "qs": q_values,
        #            "alphas": alpha_values.tolist(),
        #            "line_range": [start_line, end_line]
        #})

        start_line = end_line

with open(json_file, 'w') as jsonfile:
    json.dump(json_data, jsonfile, indent=4)