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
import gc
import faulthandler
import os
import multiprocessing

# Enable faulthandler to catch crashes
faulthandler.enable()

from ctypes_classes import CCC
from ctypes_classes import LC_irs
from analytic_classes import LightCurveAnalyzer, Complex, get_trajectory, LightCurveDTW

sys.setrecursionlimit(10000)


def run_dtw_task(a, b, th, m2, m3, source_size, lc_len, points_per_radius, q_list, alpha_list, ini_time, fin_time, result_queue):
    try:
        # Import necessary modules and classes
        from ctypes_classes import LC_irs
        from analytic_classes import LightCurveAnalyzer, Complex, get_trajectory, LightCurveDTW

        # Initialize the LC_irs object
        analyzerDTW = LightCurveDTW(a, b, th, m2, m3, source_size, lc_len, points_per_radius, "features_summary")

        line_counter = 0
        json_data = []
        for q in q_list:
          for alpha in alpha_list:
            # Run your light curve analysis
            feature_vector, light_curve, residual = analyzerDTW.run_for_params(q, alpha, ini_time, fin_time, line_counter % 23 == 0)
            line_counter += 1        
            # check if there is substantial light curve
            #if np.linalg.norm(np.array(light_curve)) > 0.5*len(light_curve):
            if sum(light_curve) > 0.5*len(light_curve):
              json_data.append({
                  "a": a,
                  "b": b,
                  "theta": th,
                  "m2": m2,
                  "m3": m3,
                  "q": q,
                  "alpha": alpha,
                  "feature_vector": feature_vector,
                  "light_curve": residual.tolist()
              })
            else:
                print("flat curve for a={a}, q={q}, alpha={alpha}")

        if len(json_data) > 0:
        # Put the result in the queue
            result_queue.put(json_data)
        else:
            result_queue.put({"error": f"No usable data for a={a}, q={q}, alpha={alpha}"})    
    

    except Exception as e:
        print(f"Process {os.getpid()} failed with exception: {e}")
        # In case of failure, send an error message or empty result
        result_queue.put({"error": f"Failed for a={a}, q={q}, alpha={alpha} with exception: {e}"})


def orchestrator_dtw():
    alphas_in_pi = 32
    #a_values = [0.8, 0.85, 0.9, 0.95, 0.98, 1.0, 1.02, 1.05, 1.07, 1.1, 1.15, 1.2]
    a_values = [0.75, 1.11, 1.12]
    b = 0.0
    q_values = [0.05, 0.1, 0.13, 0.15, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5]
    #q_values = [0.05, 0.1]
    alpha_values = np.arange(0, 2 * np.pi, np.pi / float(alphas_in_pi))
    th = 1.047197551
    m2 = 5.0e-3
    m3 = 0.0
    source_size = 1e-4
    lc_len = 512
    points_per_radius = 50
    ini_time = -1.0
    fin_time = 1.0


    # Create a queue to collect results
    result_queue = multiprocessing.Queue()

    alpha_lists = [alpha_values[i:i + alphas_in_pi] for i in range(0, len(alpha_values), min(4, alphas_in_pi))]

    json_file = "parameters_"
    start_line = 0
    end_line = 0
    
    for a in a_values:
        json_data = []
        for q in q_values:
            for alpha in alpha_values:
                p = multiprocessing.Process(target=run_dtw_task, args=(a, b, th, m2, m3, source_size, lc_len, points_per_radius, [q], [alpha], ini_time, fin_time, result_queue))
                p.start()
                p.join()  # Wait for the process to finish

                if p.exitcode != 0:
                    print(f"Task with parameters a={a}, q={q} failed with exit code {p.exitcode}. Moving to the next task.")
                else:
                    # Get the result from the queue
                    result_datum = result_queue.get()
                    
                    if "error" in result_datum:
                        print(result_datum["error"])
                    else:
                        print(f"Task with parameters a={a}, q={q}, completed successfully.")
                        json_data.extend(result_datum)
                

        try:
          # Write to JSON file after processing all q and alpha values for the current a
          with open(json_file + "a=" + str(a) + ".json", 'w') as jsonfile:
              json.dump(json_data, jsonfile, indent=2)
          print(f"JSON file created for a={a}")
        except Exception as e:
          print(f"Failed to write JSON file for a={a}: {e}")



## noisy memory leak outputs
##gc.set_debug(gc.DEBUG_LEAK)
#
## Define lens parameters.
#a = 1.01
#b = 1.0001
#theta = 1.047197551
#m2 = 5.0e-3
#m3 = 0.0
#length = 500
#
#
## Define source parameters.
#source_size = 1e-4
#lc_steps = 512
#points_per_radius = 30
#ini_time = -1.0
#fin_time = 1.0
#
#
## Define parameter ranges.
##a_values = [0.9, 0.95, 1.0, 1.1]  # Example values for a
##q_values = [0.05, 0.1, 0.15, 0.2, 0.4]  # Example values for q
##alpha_values = np.arange(0, 2 * np.pi, np.pi / 4.0)  # Example values for alpha
#
#a_values = [0.9, 0.95]  # Example values for a
#q_values = [0.05, 0.1, 0.15, 0.2, 0.3, 0.35, 0.4]  # Example values for q
##q_values = [0.05]  # Example values for q
#alpha_values = np.arange(0, 2 * np.pi, np.pi / 2.0) # Example values for alpha
##alpha_values = [np.pi / 2.0] # Example values for alpha
#
## Output files
#csv_file = "feature_matrix.csv"
#json_file = "parameters_"
#feature_file_name = "features_summary"
#
## JSON structure to hold parameter values and line ranges
#json_data = []
#
#start_line = 0
#end_line = 0


#for a in a_values:
#  analyzerDTW = LightCurveDTW(a, b, theta, m2, m3, source_size, lc_steps, points_per_radius, feature_file_name)
#
#  json_data = []
#  for q in q_values:
#      for alpha in alpha_values:
      #    try:
      #      [feature_vector, light_curve, residual] = analyzerDTW.run_for_params(q, alpha, ini_time, fin_time, end_line % 4 == 0)
      #      #writer.writerow(feature_vector)
      #      end_line += 1
      #      json_data.append({
      #          "a": a,
      #          "b": b,
      #          "theta": theta,
      #          "m2": m2,
      #          "m3": m3,
      #          "q": q,
      #          "alpha": alpha,
      #          "feature_vector": feature_vector,
      #          "light_curve": residual.tolist()
      #      })
      #    except Exception as e:
      #      print(f"\nAn error occurred while running for params: {e}")
      #      print(f"Values: a={a}, q={q}, alpha={alpha}")
      #      traceback.print_exc()


#  start_line = end_line
#
#  try:
#    # Write to JSON file after processing all q and alpha values for the current a
#    with open(json_file + "a=" + str(a) + ".json", 'w') as jsonfile:
#        json.dump(json_data, jsonfile, indent=2)
#    print(f"JSON file created for a={a}")
#  except Exception as e:
#    print(f"Failed to write JSON file for a={a}: {e}")
#
#  analyzerDTW = None
#  gc.collect()




if __name__ == "__main__":
    orchestrator_dtw()