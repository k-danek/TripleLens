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
from dtaidistance import dtw
from matplotlib.patches import Rectangle, Patch
import json

from ctypes_classes import CCC
from ctypes_classes import LC_irs

class Complex(ctypes.Structure):
   _fields_ = [("real", ctypes.c_double),
               ("imag", ctypes.c_double)
              ]

def get_single_amp_real(time, t_0, q, time_scale):
    u_r = time_scale * (time - t_0)
    u_sq = u_r * u_r + q * q
    return (u_sq + 2.0) / np.sqrt(u_sq * (u_sq + 4.0))

def get_single_amp_real_scale(time, t_0, q, amp_scale, time_scale):
    u_r = time_scale * (time - t_0)
    u_sq = u_r * u_r + q * q
    return amp_scale * (u_sq + 2.0) / np.sqrt(u_sq * (u_sq + 4.0))

def get_trajectory(q, alpha, iniTime, finTime, steps):
    trajectory = []

    for i in range (0, steps):
        time = iniTime + (finTime-iniTime)*i/(steps-1)
        point = Complex(math.cos(alpha)*time-q*math.sin(alpha), math.sin(alpha)*time+q*math.cos(alpha))
        trajectory.append(point)

    return trajectory

class LightCurveAnalyzer:
    def __init__(self, a, b, theta, m2, m3, source_size, lc_steps, points_per_radius):
        self.a = a
        self.b = b
        self.theta = theta
        self.m2 = m2
        self.m3 = m3
        self.source_size = source_size
        self.lc_steps = lc_steps
        self.points_per_radius = points_per_radius
        self.lc_irs = LC_irs(a, b, theta, m2, m3, source_size, lc_steps, points_per_radius)
        self.ccc = CCC(a, b, theta, m2, m3, 500)
        self.ccc.get_ca()
        self.cc_array = np.zeros(3000, np.cdouble)
        self.ca_array = np.zeros(3000, np.cdouble)
        self.ccc.copy_cc_ca(self.cc_array, self.ca_array)
        self.lenses_real, self.lenses_imag = self.get_lenses()

    def get_light_curve(self, q, alpha, ini_time, fin_time):
        pos_ini_x = math.cos(alpha) * ini_time - q * math.sin(alpha)
        pos_ini_y = math.sin(alpha) * ini_time + q * math.cos(alpha)
        pos_fin_x = math.cos(alpha) * fin_time - q * math.sin(alpha)
        pos_fin_y = math.sin(alpha) * fin_time + q * math.cos(alpha)

        lc_irs_array = np.zeros(self.lc_steps, np.double)

        self.lc_irs.get_lc_irs(pos_ini_x, pos_ini_y, pos_fin_x, pos_fin_y)
        self.lc_irs.copy_lc(lc_irs_array)
        
        return lc_irs_array

    def get_lenses(self):
      lens_pos = np.zeros(3, np.cdouble)
      self.ccc.copy_lenses(lens_pos)
      lenses_real = []
      lenses_imag = []
      for lens in lens_pos:
        lenses_real.append(lens.real)
        lenses_imag.append(lens.imag)
      return lenses_real, lenses_imag
    
    def get_ccc(self):
      self.ccc.get_ca()

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
      
      return cc_real, cc_imag, ca_real, ca_imag


    def fit_single_light_curve(self, time_series, initial_guesses):
        peaks, _ = find_peaks(time_series, width=(2, int(self.lc_steps / 5)))
        time_series_no_peaks = np.copy(time_series)
        time_series_no_peaks[peaks] = np.nan  # Remove peaks
        mask = ~np.isnan(time_series_no_peaks)
        x_data = np.arange(len(time_series))
        
        # Initial fit
        popt1, pcov1 = curve_fit(get_single_amp_real, x_data[mask], time_series_no_peaks[mask], p0=initial_guesses)
        fitted_smooth1 = get_single_amp_real(x_data, *popt1)
        
        # Calculate residuals and standard deviation
        residuals1 = time_series - fitted_smooth1
        std_devs = np.abs(residuals1)
        
        # Mask the 10% of points with the highest standard deviation
        threshold = np.percentile(std_devs, 90)
        mask_high_std = std_devs < threshold
        
        # Second fit
        popt2, pcov2 = curve_fit(get_single_amp_real, x_data[mask & mask_high_std], time_series_no_peaks[mask & mask_high_std], p0=initial_guesses)
        fitted_smooth2 = get_single_amp_real(x_data, *popt2)
        
        residuals2 = time_series - fitted_smooth2
        
        return (x_data, time_series, fitted_smooth1, fitted_smooth2, residuals2, popt1, popt2, peaks)

    def plot_results(self, x_data, time_series, fitted_smooth1, fitted_smooth2, residuals2, peaks, popt1, popt2, q, alpha, ini_time, fin_time, filename):
        fit_label1 = f'Initial fit values:\nt_0={popt1[0]:.2f}, q={popt1[1]:.2f}, time_scale_reversed={1.0/popt1[2]:.2f}'
        fit_label2 = f'Second fit values:\nt_0={popt2[0]:.2f}, q={popt2[1]:.2f}, time_scale_reversed={1.0/popt2[2]:.2f}'
        
        # Generate trajectory
        trajectory = get_trajectory(q, alpha, ini_time, fin_time, self.lc_steps)
        traj_real = [p.real for p in trajectory]
        traj_imag = [p.imag for p in trajectory]
        
        fig, ax = plt.subplots(2, 2, figsize=(12, 8), gridspec_kw={'height_ratios': [3, 1]})
        
        # Main plot
        ax[0, 0].plot(x_data, time_series, label='Original Time Series')
        #ax[0, 0].plot(x_data, fitted_smooth1, label='Initial Fitted Smooth Background')
        ax[0, 0].plot(x_data, fitted_smooth2, label='Second Fitted Smooth Background')
        ax[0, 0].scatter(x_data[peaks], time_series[peaks], color='red', label='Detected Peaks')
        ax[0, 0].legend()
        ax[0, 0].set_xlabel('Time')
        ax[0, 0].set_ylabel('Value')
        ax[0, 0].set_title(f'Time Series with Fitted Smooth Backgrounds\n{fit_label1}\n{fit_label2}')
        
        # Residuals plot
        ax[1, 0].plot(x_data, residuals2, label='Residuals (Second Fit)')
        ax[1, 0].set_xlabel('Time')
        ax[1, 0].set_ylabel('Residuals')
        ax[1, 0].legend()
        
        # Critical curve and caustic plot
        ax[0, 1].scatter(self.cc_array.real, self.cc_array.imag, color='blue', s=1, label='Critical Curve')
        ax[0, 1].scatter(self.ca_array.real, self.ca_array.imag, color='red', s=1, label='Caustic')
        ax[0, 1].scatter(self.lenses_real, self.lenses_imag, color='green', s=20, label='Lenses')
        ax[0, 1].plot(traj_real, traj_imag, color='black', label='Trajectory')
        ax[0, 1].legend()
        ax[0, 1].set_xlabel('Real Part')
        ax[0, 1].set_ylabel('Imaginary Part')
        ax[0, 1].set_title('Critical Curve and Caustic')
        
        plt.tight_layout()
        plt.savefig(filename+".png")
        plt.close()
        
        print("Plot saved as 'time_series_analysis_with_residuals.png'")
        print(f"Initial fit parameters: t_0={popt1[0]:.2f}, q={popt1[1]:.2f}, time_scale_reversed={1.0/popt1[2]:.2f}")
        print(f"Second fit parameters: t_0={popt2[0]:.2f}, q={popt2[1]:.2f}, time_scale_reversed={1.0/popt2[2]:.2f}")
    
    def print_residuals(self, x_data, residuals2, filename):
      with open(filename+".txt", 'w') as f:
        #f.write("Residuals (Second Fit):\n")
          for residual in residuals2:
            f.write(f"{residual}\n")


class Feature:
  def __init__(self, feature_name, feature_data, feature_group):
    self.feature_name = feature_name
    self.feature_data = feature_data
    self.feature_group = feature_group

class Feature_group:
  def __init__(self, feature_group_name, feature_names, feature_overlap_group,  group_index):
    self.feature_group_name = feature_group_name
    self.feature_names = feature_names
    self.feature_overlap_group = feature_overlap_group
    self.group_index = group_index


class LightCurveDTW:
    def __init__(self, a, b, theta, m2, m3, source_size, lc_steps, points_per_radius, feature_file_name):
        self.analyzer = LightCurveAnalyzer(a, b, theta, m2, m3, source_size, lc_steps, points_per_radius)

        self.lc_steps = lc_steps
        self.points_per_radius = points_per_radius
        self.feature_names = [] # array of string feature names
        self.lc_irs = LC_irs(a, b, theta, m2, m3, source_size, lc_steps, points_per_radius)
        self.short_signals = {}
        self.feature_groups = {}
        with open(f"./features/{feature_file_name}.json", 'r') as f:
            feature_data = json.load(f)
        
            for feature_group in feature_data:
              feature_group_name = feature_group["group_name"]
              feature_group_names = feature_group["names"]
              feature_ovelap_group = feature_group["related_groups"]
              feature_group_index = feature_group["group_index"]
              self.feature_groups[feature_group_name] = Feature_group(feature_group_name, feature_group_names, feature_ovelap_group, feature_group_index)
              for feature_name in feature_group_names:
                short_feature_data = np.loadtxt(f"./features/{feature_name}.txt")
                self.short_signals[feature_name] = Feature(feature_name, short_feature_data, self.feature_groups[feature_group_name])


        print(self.short_signals)

        # Color map for different short signals
        #cmap = plt.get_cmap('tab9')
        cmap = plt.get_cmap('Accent')
        unique_feature_groups = list(self.feature_groups.keys())
        #self.colors = [cmap(i) for i in range(len(feature_data))]
        self.color_map = {group: cmap(i / len(unique_feature_groups)) for i, group in enumerate(unique_feature_groups)}

    def get_color_for_group(self, feature_group_name):
        # Retrieve the color for a specific feature group
        return self.color_map.get(feature_group_name, (0, 0, 0, 1))

    def sliding_window_dtw(self, short_signal, long_signal, window_size):
        distances = []
        positions = []
        alignments = []
        step_size = max(1, window_size // 4)
        short_signal_norm = np.linalg.norm(short_signal)

        for i in range(0, len(long_signal) - window_size + 1, step_size):
            window = long_signal[i:i + window_size]
            window_norm = np.linalg.norm(window)

            # Only calculate DTW distance if window has significant norm
            if window_norm >= 0.4 * short_signal_norm:
                normalized_short_signal = short_signal / short_signal_norm
                normalized_window = window / window_norm
                dist, paths = dtw.warping_paths(normalized_short_signal, normalized_window)
                best_path = dtw.best_path(paths)
                distances.append(dist)
                positions.append(i)
                alignments.append((normalized_short_signal * window_norm, normalized_window, best_path))

        return distances, positions, alignments


    #def get_cleaned_overlapping_features_II(self, matched_tuples):
    # checks if two fitted features are overlapping
    def is_overlapping(self, matched_tuple_1, matched_tuple_2):
      pos_1, len_1 =  matched_tuple_1[1], len(matched_tuple_1[2].feature_data)-1
      pos_2, len_2 =  matched_tuple_2[1], len(matched_tuple_2[2].feature_data)-1
      type_1, overlap_types_1 = matched_tuple_1[2].feature_group.feature_group_name, matched_tuple_1[2].feature_group.feature_overlap_group 
      type_2, overlap_types_2 = matched_tuple_2[2].feature_group.feature_group_name, matched_tuple_2[2].feature_group.feature_overlap_group 

      do_overlap = False

      if (pos_1 <= pos_2+len_2-1) and (pos_2 <= pos_1+len_1-1):
        if (type_1 in overlap_types_2) and (type_2 in overlap_types_1) :
          do_overlap = True
      
      return do_overlap


    def is_bested_by(self, idx_1, matched_tuples):
      do_best_feature = set()
      for idx_2 in range(0,len(matched_tuples)):
        if idx_2 != idx_1:
          if self.is_overlapping(matched_tuples[idx_1], matched_tuples[idx_2]):
            if matched_tuples[idx_1][0] > matched_tuples[idx_2][0]:
              do_best_feature.add(idx_1)

      return do_best_feature

    # removes features that overlap with features that we for sure want to keep
    def remove_overlap_layer(self, overlap_list, cleaned_set, drop_set):
      layer_removed = False
      
      # collect cleaned set, then collect drop se, i.e., the set which cleaned set overlaps
      for idx_1 in range(0,len(overlap_list)):
        if len(overlap_list[idx_1]) == 0:
          # the tuples that are not bested, can be cleaned
          cleaned_set.add(idx_1)
          # the ones that are cleaned, will remove those that are not
          for idx_2 in range(0,len(overlap_list)):
            if idx_2 != idx_1:
              if idx_1 in overlap_list[idx_2]:
                drop_set.add(idx_2)

      # clean those overlaps that are in drop set      
      for idx_1 in range(0,len(overlap_list)):
        if idx_1 not in cleaned_set:
          for idx_2 in drop_set:
            if idx_2 in overlap_sets[idx_1]:
              overlap_list[idx_1].remove(idx_2)
              layer_removed = True

      return overlap_list, cleaned_set, drop_set, layer_removed

    def get_cleaned_overlapping_features_II(self, matched_tuples):
      drop_set = set()
      cleaned_set = set()

      # Sort by position to make overlaps less confusing
      matched_tuples = sorted(matched_tuples, key=lambda x: x[1])

      overlap_list = []

      # collect sets
      for idx_1 in range(0,len(matched_tuples)):
        overlap_list.append(self.is_bested_by(idx_1, matched_tuples))


      should_remove_layer = True
      while should_remove_layer:
        overlap_list, cleaned_set, drop_set, should_remove_layer  = self.remove_overlap_layer(overlap_list, cleaned_set, drop_set)
      
      return [matched_tuples[i] for i in cleaned_set if i not in drop_set]


    # Takes array of matches and removes those that have too high threshold
    def clean_distant_features(self, matched_tuples, threshold):
      close_enough_matches = [] 
      for match in matched_tuples:
        if match[0] < threshold:
          close_enough_matches.append(match)

      return close_enough_matches 

    def run_for_params(self, q, alpha, ini_time, fin_time, plot_fig):
        time_series = self.analyzer.get_light_curve(q, alpha, ini_time, fin_time)
        initial_guesses = [self.lc_steps/2.0, q, 1.0/float(self.lc_steps)]
        
        # Generate trajectory
        trajectory = get_trajectory(q, alpha, ini_time, fin_time, self.lc_steps)
        traj_real = [p.real for p in trajectory]
        traj_imag = [p.imag for p in trajectory]
        
        x_data, time_series, fitted_smooth1, fitted_smooth2, residuals, popt1, popt2, peaks = self.analyzer.fit_single_light_curve(time_series, initial_guesses)
        
        best_matches = []
        
        for feature_name, short_signal in self.short_signals.items():
            window_size = len(short_signal.feature_data)
            distances, positions, alignments = self.sliding_window_dtw(short_signal.feature_data, residuals, window_size)
            min_distances_indices = np.argsort(distances)[:6]  # Get top 4 matches for each feature

            feature_color = self.get_color_for_group(short_signal.feature_group)

            for idc in min_distances_indices:
              best_matches.append((distances[idc], positions[idc], short_signal, feature_color, alignments[idc]))

        print("\n\n len of best matches before distance cleaning "+str(len(best_matches)))

        best_matches = self.clean_distant_features(best_matches, 0.7)

        print("\n\n len of best matches before overlap cleaning "+str(len(best_matches)))
        # In order not to introduce more features than there already is, I took the best distance from overlapping features
        best_matches = self.get_cleaned_overlapping_features_II(best_matches)


        # Sort all matches by distance and take only the first 12

        print("\n\n len of best matches before sorting position cleaning "+str(len(best_matches)))
        
        # Sort by position
        best_matches = sorted(best_matches, key=lambda x: x[1])
        
        #########################################################
        if plot_fig:
          try:
            # Calculate the number of rows needed for short signals
            n_short_signal_rows = len(best_matches)
            n_main_rows = 2  # Two main rows for time series and residuals

            fig = plt.figure(figsize=(24, 12))  # Adjusted figure size to accommodate new plot
            gs = fig.add_gridspec(n_main_rows, 3, width_ratios=[2, 3, 1])

            # Plot the additional rectangular plot
            ax0 = fig.add_subplot(gs[:, 0])
            ax0.scatter(self.analyzer.cc_array.real, self.analyzer.cc_array.imag, color='blue', s=1, label='Critical Curve')
            ax0.scatter(self.analyzer.ca_array.real, self.analyzer.ca_array.imag, color='red', s=1, label='Caustic')
            ax0.scatter(self.analyzer.lenses_real, self.analyzer.lenses_imag, color='green', s=20, label='Lenses')
            ax0.plot(traj_real, traj_imag, color='black', label='Trajectory')
            ax0.legend()
            ax0.set_xlabel('Real Part')
            ax0.set_ylabel('Imaginary Part')
            ax0.set_title('Critical Curve and Caustic')

            # Plot time series
            ax1 = fig.add_subplot(gs[0, 1])
            ax1.plot(time_series, label='Time Series')
            ax1.set_title(f'Time Series for q={q}, alpha={alpha}')
            ax1.legend()

            legend_patches = []  # Store patches for legend

            # Plot residuals
            ax2 = fig.add_subplot(gs[1, 1])
            ax2.plot(residuals, label='Residuals')
            for dist, pos, short_signal, color, (normalized_short_signal, normalized_window, best_path) in best_matches:
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
                legend_patches.append(Patch(facecolor=color, edgecolor=color, label=f'{short_signal.feature_name}; d: {dist:.3f}; p: {pos:d}'))

            ax2.set_title('Residuals with Matches')
            ax2.legend()
            ax2.legend(handles=legend_patches, loc='upper right')

            # Add new grid spec for short signals
            gs_short = fig.add_gridspec(n_short_signal_rows, 1, left=0.85, right=0.99, hspace=0.4)

            #:# Plot short signals
            #for i, short_signal in enumerate(self.short_signals.values()):
            #    ax = fig.add_subplot(gs_short[i, 0])
            #    ax.plot(short_signal.feature_data, color=self.colors[i % len(self.colors)])
            #    ax.set_title(short_signal.feature_name)

            plt.tight_layout()
            plt.savefig(f"q_{q}_alpha_{alpha}.png")
            plt.close(fig)
          except Exception as e:
            print(f"\nAn error occurred while creating or saving the plot: {e}")
          #finally:
          #  plt.close(fig)
        #######################################################################        


        feature_vector = []
        for match in best_matches:
          feature_vector.append(float(match[1]))
          feature_vector.append(match[2].feature_group.group_index)
    
        while len(feature_vector) < 15:
          feature_vector.append(-1)
          feature_vector.append(-1.0)


        return [feature_vector, time_series, residuals] 