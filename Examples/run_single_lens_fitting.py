import ctypes
from ctypes import *
import time
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import math
import sys

from ctypes_classes import CCC
from ctypes_classes import LC_irs

from analytic_classes import LightCurveAnalyzer,  get_single_amp_real, get_single_amp_real_scale, get_trajectory

sys.setrecursionlimit(500000) 

# Parameters

# Define lens parameters.
a = 1.01
b = 1.0001
theta = 1.047197551
m2 = 5.0e-3
#m3 = 1e-5
m3 = 0.0
length = 500


# Define source parameters.
source_size = 1e-4
lc_steps = 512
points_per_radius = 30
ini_time = -1.0
fin_time = 1.0

q_values = [0.05, 0.1, 0.12, 0.15, 0.18, 0.2, 0.25, 0.3]  # Example values for q
#alpha_values = [np.pi / 8.0, np.pi / 5.0, np.pi / 2.0, 5.0 * np.pi / 8.0, 3.0 * np.pi / 4.0, 7.0 * np.pi / 8.0]  # Example values for alpha
alpha_values =  np.arange(0, 2 * np.pi, np.pi / 8.0)  # Example values for alpha

analyzer = LightCurveAnalyzer(a, b, theta, m2, m3, source_size, lc_steps, points_per_radius)

for q in q_values:
    for alpha in alpha_values:
        time_series = analyzer.get_light_curve(q, alpha, ini_time, fin_time)
        initial_guesses = [lc_steps/2.0, q, 1.0/float(lc_steps)]
        x_data, time_series, fitted_smooth1, fitted_smooth2, residuals2, popt1, popt2, peaks = analyzer.fit_single_light_curve(time_series, initial_guesses)
        filename = f"single_lens_fit_q={q:.2f}_a={alpha/np.pi:.2f}"
        analyzer.plot_results(x_data, time_series, fitted_smooth1, fitted_smooth2, residuals2, peaks, popt1, popt2, q, alpha, ini_time, fin_time, filename)
        analyzer.print_residuals(x_data, residuals2, filename)
