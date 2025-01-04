# This script takes a range of parameters, generate light curves, labels with feature vectors

from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import sys
import json
import faulthandler
import os
import multiprocessing

# Enable faulthandler to catch crashes
faulthandler.enable()

sys.setrecursionlimit(10000)


def run_dtw_task(a, b, th, m2, m3, source_size, lc_len, points_per_radius, q_list, alpha_list, ini_time, fin_time, result_queue):
    try:
        # Import necessary modules and classes
        from examples.analytic_classes import LightCurveDTW

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
            result_queue.put({"error": f"No usable data for a={a}"})    
    

    except Exception as e:
        print(f"Process {os.getpid()} failed with exception: {e}")
        # In case of failure, send an error message or empty result
        result_queue.put({"error": f"Failed for a={a} with exception: {e}"})


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

    #alpha_lists = [alpha_values[i:i + alphas_in_pi] for i in range(0, len(alpha_values), min(4, alphas_in_pi))]

    json_file = "./ml_project/parameters_"
    
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

if __name__ == "__main__":
    orchestrator_dtw()