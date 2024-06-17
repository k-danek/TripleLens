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
from run_single_lens_fitting import LightCurveAnalyzer

sys.setrecursionlimit(500000) 

def sliding_window_dtw(short_signal, long_signal, window_size):
    distances = []
    positions = []
    alignments = []
    step_size = max(1, window_size // 2)
    short_signal_norm = np.linalg.norm(short_signal)
    
    for i in range(0, len(long_signal) - window_size + 1, step_size):
        window = long_signal[i:i + window_size]
        window_norm = np.linalg.norm(window)
        
        # Only calculate DTW distance if window has significant norm
        if window_norm >= 0.1 * short_signal_norm:
            normalized_short_signal = short_signal / short_signal_norm
            normalized_window = window / window_norm
            dist, paths = dtw.warping_paths(normalized_short_signal, normalized_window)
            best_path = dtw.best_path(paths)
            distances.append(dist)
            positions.append(i)
            alignments.append((normalized_short_signal * window_norm, normalized_window, best_path))
    
    return distances, positions, alignments

def group_overlapping_windows(distances, positions, window_size):
    sorted_indices = np.argsort(positions)
    grouped_windows = []
    current_group = []

    for idx in sorted_indices:
        pos = positions[idx]
        if not current_group or pos <= current_group[-1][1]:
            current_group.append((distances[idx], pos, pos + window_size - 1))
        else:
            grouped_windows.append(current_group)
            current_group = [(distances[idx], pos, pos + window_size - 1)]
    
    if current_group:
        grouped_windows.append(current_group)
    
    return grouped_windows

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

analyzer = LightCurveAnalyzer(a, b, theta, m2, m3, source_size, lc_steps, points_per_radius)

#for q in q_values:
#    for alpha in alpha_values:
#        time_series = analyzer.get_light_curve(q, alpha, ini_time, fin_time) 
#        initial_guesses = [lc_steps/2.0, q, 1.0/float(lc_steps)]
#        x_data, time_series, fitted_smooth1, fitted_smooth2, residuals, popt1, popt2, peaks = analyzer.fit_light_curve(time_series, initial_guesses)

feature_names = ["caustic_entry", "caustic_exit", "cusp_approach_a", "cusp_approach_b", "cusp_transversal", "dip", "double_crossing"]
#feature_names = ["caustic_entry", "caustic_exit", "cusp_transversal", "dip", "double_crossing"]
short_signals = {}
for feature in feature_names:
    short_signals[feature] = np.loadtxt(f"./features/{feature}.txt")

# Color map for different short signals
#colors = plt.cm.get_cmap('tab10', len(feature_names))
#colors = plt.cm.tab10.colors

# Color map for different short signals
cmap = plt.get_cmap('tab10')
colors = [cmap(i) for i in range(len(feature_names))]

for q in q_values:
    for alpha in alpha_values:
        # Generate the long signal
        time_series = analyzer.get_light_curve(q, alpha, ini_time, fin_time)
        initial_guesses = [lc_steps/2.0, q, 1.0/float(lc_steps)]
        x_data, time_series, fitted_smooth1, fitted_smooth2, residuals, popt1, popt2, peaks = analyzer.fit_light_curve(time_series, initial_guesses)
        
        best_matches = []
        
        for idx, feature_name in enumerate(feature_names):
            short_signal = short_signals[feature_name]
            window_size = len(short_signal)
            distances, positions, alignments = sliding_window_dtw(short_signal, residuals, window_size)
            min_distances_indices = np.argsort(distances)[:4]  # Get top 4 matches
            for idc in min_distances_indices:
                best_matches.append((distances[idc], positions[idc], feature_name, colors[idx], alignments[idc]))

        # Sort all matches by distance
        best_matches = sorted(best_matches)[:5]

        # Calculate the number of rows needed for short signals
        n_short_signal_rows = len(feature_names)
        n_main_rows = 2  # Two main rows for time series and residuals

        fig = plt.figure(figsize=(18, 12))
        gs = fig.add_gridspec(n_main_rows, 2, width_ratios=[3, 1])

        # Plot time series
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(time_series, label='Time Series')
        ax1.set_title(f'Time Series for q={q}, alpha={alpha}')
        ax1.legend()

        legend_patches = []  # Store patches for legend

        # Plot residuals
        ax2 = fig.add_subplot(gs[1, 0])
        ax2.plot(residuals, label='Residuals')
        for dist, pos, feature_name, color, (normalized_short_signal, normalized_window, best_path) in best_matches:
            #ax2.axvline(x=pos, linestyle='--', color=color, label=f'{feature_name}; d: {dist:.3f}')
            aligned_signal = np.zeros_like(residuals)
            for (i, j) in best_path:
                aligned_signal[pos + j] = normalized_short_signal[i]
            ax2.plot(aligned_signal, color=color)
            # Calculate bounding box
            aligned_indices = [pos + j for (i, j) in best_path]
            min_idx = min(aligned_indices)
            max_idx = max(aligned_indices)
            min_val = min(aligned_signal[aligned_indices])
            max_val = max(aligned_signal[aligned_indices])
            rect = Rectangle((min_idx, min_val), max_idx - min_idx, max_val - min_val, linewidth=2, edgecolor=color, facecolor='none')
            ax2.add_patch(rect)
            legend_patches.append(Patch(facecolor=color, edgecolor=color, label=f'{feature_name}; d: {dist:.3f}'))

        ax2.set_title('Residuals with Matches')
        ax2.legend()
        ax2.legend(handles=legend_patches, loc='upper right')

        # Add new grid spec for short signals
        gs_short = fig.add_gridspec(n_short_signal_rows, 1, left=0.75, right=0.95, hspace=0.4)

        # Plot short signals
        for i, feature_name in enumerate(feature_names):
            ax = fig.add_subplot(gs_short[i, 0])
            ax.plot(short_signals[feature_name], color=colors[i % len(colors)])
            ax.set_title(feature_name)

        plt.tight_layout()
        plt.savefig(f"q_{q}_alpha_{alpha}.png")
        plt.close(fig)