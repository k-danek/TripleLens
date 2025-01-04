# This stript picks already pre-clustered JSON data and prints CCC with coresponding residuals
# It serves to observe any regularity in the clusters. 
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import sys
import multiprocessing
import queue

from examples.analytic_classes import get_trajectory

sys.setrecursionlimit(1500000)


def get_ccc_orchestrated(a, b, th, m2, m3, q, alpha, lc_len, source_size, points_per_radius):
    # Create a queue to collect results
    manager = multiprocessing.Manager()
    result_queue = manager.Queue()
    p = multiprocessing.Process(target=run_ccc_task, args=(a, b, th, m2, m3, source_size, lc_len, points_per_radius, q, alpha, result_queue))
    timeout_value = 5

    p.start()

    try:
        # Wait for the result with a timeout
        result_datum = result_queue.get(timeout=timeout_value)
    except queue.Empty:
        print(f"Timeout: Task with parameters a={a}, q={q} did not return a result in time.")
        p.terminate()
        p.join()
        return False, None

    p.join(timeout=timeout_value)

    if p.is_alive():
        print(f"Task with parameters a={a}, q={q} did not finish in time. Terminating.")
        p.terminate()
        p.join()

    if p.exitcode != 0:
        print(f"Task with parameters a={a}, q={q} failed with exit code {p.exitcode}. Moving to the next task.")
        return False, None

    if "error" in result_datum:
        print(result_datum["error"])
        return False, None
    else:
        print(f"Task with parameters a={a}, q={q}, completed successfully.")
        return True, result_datum

def run_ccc_task(a, b, th, m2, m3, source_size, lc_len, points_per_radius, q, alpha, result_queue):
    try:
        # Import necessary modules and classes
        from examples.analytic_classes import LightCurveAnalyzer 

        analyzer = LightCurveAnalyzer(a, b, th, m2, m3, source_size, lc_len, points_per_radius)

        if len(analyzer.cc_array) > 0:
        # Put the result in the queue
          result_queue.put({
            "cc_real": analyzer.cc_array.real,
            "cc_imag": analyzer.cc_array.imag,
            "ca_real": analyzer.ca_array.real,
            "ca_imag": analyzer.ca_array.imag,
            "lenses_real": analyzer.lenses_real, 
            "lenses_imag": analyzer.lenses_imag 
          })

        else:
          result_queue.put({"error": f"No usable data for a={a}, q={q}, alpha={alpha}"}) 

    except Exception as e:
        print(f"Process {os.getpid()} failed with exception: {e}")
        # In case of failure, send an error message or empty result
        result_queue.put({"error": f"Failed for a={a}, q={q}, alpha={alpha} with exception: {e}"})

def plot_event_cc(a, b, theta, m2, m3, q, alpha, feature_vector, light_curve, subdirectory):
    lc_steps = 512
    points_per_radius = 30

    ini_time = -1.0
    fin_time = 1.0

    x_data = range(lc_steps)

    succeeded, result = get_ccc_orchestrated(a, b, theta, m2, m3, q, alpha, lc_steps, 1.0e-4, points_per_radius)
    if not succeeded:
      succeeded, result = get_ccc_orchestrated(a, b, theta, m2, m3, q, alpha, lc_steps, 1.0e-4, points_per_radius)
      if not succeeded:
        print("Failed on second try")
        return

    # Generate trajectory
    trajectory = get_trajectory(q, alpha, ini_time, fin_time, lc_steps)
    traj_real = [p.real for p in trajectory]
    traj_imag = [p.imag for p in trajectory]
    
    _, ax = plt.subplots(2, 1, figsize=(4, 8), gridspec_kw={'height_ratios': [2, 1]})

    # Residuals plot
    residual = ax[1].plot(x_data, light_curve, label='Residuals (Second Fit)')
    ax[1].set_xlabel('Time')
    ax[1].set_ylabel('Residuals')
    ax[1].legend()

    # Initialize handles and labels with the residuals plot
    handles = [residual]
    labels = ['Residuals (Second Fit)']
    colors = ['r', 'g', 'b', 'y', 'c', 'm']

    # Add vertical lines for the feature vector
    for i in range(3, len(feature_vector)-2, 3):
        position = feature_vector[i]
        label_value = feature_vector[i+1]
        if label_value != -1:
            line = ax[1].axvline(x=position, color=colors[int(i/3)], linestyle='--', label=f"Feature {int(label_value)}")
            handles.append(line)
            labels.append(f"Feature {int(label_value)}")

    # Set the legend with all handles and labels
    ax[1].legend(handles, labels)
    ax[1].set_title(f"t_0={feature_vector[0]:.1f}, q={feature_vector[1]:.1e}, t_scale={feature_vector[2]:.1e}")

    # Critical curve and caustic plot
    ax[0].scatter(result["cc_real"], result["cc_imag"], color='blue', s=1, label='Critical Curve')
    ax[0].scatter(result["ca_real"], result["ca_imag"], color='red', s=1, label='Caustic')
    ax[0].scatter(result["lenses_real"], result["lenses_imag"], color='green', s=20, label='Lenses')
    ax[0].plot(traj_real, traj_imag, color='black', label='Trajectory')
    ax[0].legend()
    ax[0].set_xlabel('Real Part')
    ax[0].set_ylabel('Imaginary Part')
    ax[0].set_title('Critical Curve and Caustic')


    filename = os.path.join(subdirectory, f"cluster_{cluster_num:d}_a={a:.3f}_q={q:.2f}_alpha={alpha/np.pi:.4f}.png")

    os.makedirs(subdirectory, exist_ok=True)

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    

def plot_cluster_file(file_path, subdirectory):
    with open(file_path, 'r') as f:
        data = json.load(f) 
    for event in data:
        a = event["a"] 
        b = event["b"] 
        theta = event["theta"]
        m2 = event["m2"]
        m3 = event["m3"]
        q = event["q"]
        alpha = event["alpha"]
        feature_vector = event["feature_vector"]
        light_curve = event["light_curve"]
        plot_event_cc(a, b, theta, m2, m3, q, alpha, feature_vector, light_curve, subdirectory) 


if __name__ == "__main__":
  num_clusters = [2,3,4,5,6]
  num_of_subcluster_dict = {2: 4, 3: 12, 4: 6, 5: 1, 6: 1}

  output_dir = os.path.join('ml_project', 'cluster_results')
  os.makedirs(output_dir, exist_ok=True)
  
  for cluster_num in num_clusters:
    if cluster_num in num_of_subcluster_dict:  
      for subcluster_num in range(num_of_subcluster_dict[cluster_num]):
        print("\n\n\n**********cluster_num:" + str(cluster_num)+"**********\n\n\n")
        subdirectory = os.path.join(output_dir, f"{cluster_num}_{subcluster_num}")
        file_path = os.path.join(output_dir, f'subcluster_{cluster_num}_{subcluster_num}.json')
        print("file path:"+str(file_path))
        plot_cluster_file(file_path, subdirectory)

