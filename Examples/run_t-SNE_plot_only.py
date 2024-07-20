import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import os
import json
import sys

from analytic_classes import LightCurveAnalyzer,  get_single_amp_real, get_single_amp_real_scale, get_trajectory

sys.setrecursionlimit(1500000)


def plot_event(a, b, theta, m2, m3, q, alpha):
    lc_steps = 512
    points_per_radius = 30
    ini_time = -1.0
    fin_time = 1.0
    analyzer = LightCurveAnalyzer(a, b, theta, m2, m3, 1.0e-4, lc_steps, points_per_radius)
    
    try:
        time_series = analyzer.get_light_curve(q, alpha, ini_time, fin_time)
    except Exception as e:
        print(f"An error occurred while processing file {file_path}: {e}")
    
    initial_guesses = [lc_steps/2.0, q, 1.0/float(lc_steps)]
    x_data, time_series, fitted_smooth1, fitted_smooth2, residuals2, popt1, popt2, peaks = analyzer.fit_single_light_curve(time_series, initial_guesses)
    filename = f"cluster_{cluster_num:d}_fit_a={a:.3f}_q={q:.2f}_a={alpha/np.pi:.4f}"
    analyzer.plot_results(x_data, time_series, fitted_smooth1, fitted_smooth2, residuals2, peaks, popt1, popt2, q, alpha, ini_time, fin_time, filename)
    del analyzer
    del time_series

def plot_event_precalculated(a, b, theta, m2, m3, q, alpha, cluster_name, feature_vector, light_curve):
    # Create the title for the plot
    title = f"a={a}, b={b}, theta={theta}, m2={m2:.3e}, m3={m3}, q={q:.4f}, alpha={alpha:.4f}"

    # Create the filename
    cluster_name = os.path.basename(cluster_name).rsplit('.', 1)[0]
    filename = f"{cluster_name}_a{a:.3f}_m2{m2:.3e}_q{q:.3f}_alpha{alpha:.3f}.png"

    # Plot the light curve
    plt.figure(figsize=(10, 6))
    plt.plot(light_curve, label="Light Curve")

    # Add vertical lines for the feature vector
    for i in range(0, len(feature_vector), 2):
        position = feature_vector[i]
        label = feature_vector[i+1]
        if label != -1:
            plt.axvline(x=position, color='r', linestyle='--', label=f"Feature {label}")

    # Add title and labels
    plt.title(title)
    plt.xlabel("Time")
    plt.ylabel("Intensity")
    plt.legend()

    # Save the plot as a PNG file
    plt.savefig(filename)
    plt.close()


def plot_cluster_file(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f) 
    for event in data:
        #print("event:"+str(event))
        a = event["a"] 
        b = event["b"] 
        theta = event["theta"]
        m2 = event["m2"]
        m3 = event["m3"]
        q = event["q"]
        alpha = event["alpha"]
        feature_vector = event["feature_vector"]
        light_curve = event["light_curve"]
        plot_event_precalculated(a, b, theta, m2, m3, q, alpha, file_path, feature_vector, light_curve) 


num_clusters = [0,1,2,3,4,5]

output_dir = 'cluster_results'

for cluster_num in num_clusters:
    print("\n\n\n**********cluster_num:" + str(cluster_num)+"**********\n\n\n")
    file_path = os.path.join(output_dir, f'cluster_{cluster_num}.json')
    print("file path:"+str(file_path))
    plot_cluster_file(file_path)

