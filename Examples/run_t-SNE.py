import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import os
import json

from analytic_classes import LightCurveAnalyzer,  get_single_amp_real, get_single_amp_real_scale, get_trajectory

def plot_event(a, b, theta, m2, m3, q, alpha):
    lc_steps = 512
    points_per_radius = 30
    ini_time = -1.0
    fin_time = 1.0
    analyzer = LightCurveAnalyzer(a, b, theta, m2, m3, 1.0e-4, lc_steps, points_per_radius)
    time_series = analyzer.get_light_curve(q, alpha, ini_time, fin_time)
    initial_guesses = [lc_steps/2.0, q, 1.0/float(lc_steps)]
    x_data, time_series, fitted_smooth1, fitted_smooth2, residuals2, popt1, popt2, peaks = analyzer.fit_single_light_curve(time_series, initial_guesses)
    filename = f"cluster_{cluster_num:d}_fit_a={a:.3f}_q={q:.2f}_a={alpha/np.pi:.4f}"
    analyzer.plot_results(x_data, time_series, fitted_smooth1, fitted_smooth2, residuals2, peaks, popt1, popt2, q, alpha, ini_time, fin_time, filename)

def plot_cluster_file(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f) 
    for event in data:
        print("event:"+str(event))
        a = event["a"] 
        b = event["b"] 
        theta = event["theta"]
        m2 = event["m2"]
        m3 = event["m3"]
        q = event["q"]
        alpha = event["alpha"]
        plot_event(a, b, theta, m2, m3, q, alpha) 



# Read the JSON file
file_path = 'parameters_all_768.json'
with open(file_path, 'r') as f:
    data = json.load(f)

# Extract feature vectors from the JSON structure
feature_vectors = [item['feature_vector'] for item in data]

# Convert to DataFrame for easier handling
df = pd.DataFrame(feature_vectors)

# Fill NaN values with the column mean
df.fillna(df.mean(), inplace=True)

# Extract feature matrix
feature_matrix = df.values

# Perform t-SNE
tsne = TSNE(n_components=2, random_state=42)
tsne_results = tsne.fit_transform(feature_matrix)

# Perform k-means clustering on the t-SNE results
num_clusters = 6  # You can change this number based on your needs
kmeans = KMeans(n_clusters=num_clusters, random_state=42)
clusters = kmeans.fit_predict(tsne_results)

# Create a directory for cluster output files
output_dir = 'cluster_results'
os.makedirs(output_dir, exist_ok=True)

lc_steps = 512
points_per_radius = 30
ini_time = -1.0
fin_time = 1.0

# Output the cluster results to separate files
for cluster_num in range(num_clusters):
    cluster_indices = np.where(clusters == cluster_num)[0]
    cluster_data = [data[idx] for idx in cluster_indices]
    with open(os.path.join(output_dir, f'cluster_{cluster_num}.json'), 'w') as f:
        cluster_data = [data[idx] for idx in cluster_indices]
        json.dump(cluster_data, f, indent=4)

print("\n\n\n********* outputted clusters *********\n\n\n")


#for cluster_num in range(num_clusters):
#    print("\n\n\n**********cluster_num:" + str(cluster_num)+"**********\n\n\n")
#    file_path = os.path.join(output_dir, f'cluster_{cluster_num}.json')
#    print("file path:"+str(file_path))
    #plot_cluster_file(file_path)


# Plot the t-SNE results with clusters
plt.figure(figsize=(10, 6))
colors = plt.cm.get_cmap('tab10', num_clusters)
for cluster_num in range(num_clusters):
    indices = np.where(clusters == cluster_num)[0]
    plt.scatter(tsne_results[indices, 0], tsne_results[indices, 1], 
                c=[colors(cluster_num)], s=50, alpha=0.7, label=f'Cluster {cluster_num}')
    

    
plt.title('t-SNE Visualization of Feature Vectors with Clustering')
plt.xlabel('t-SNE Component 1')
plt.ylabel('t-SNE Component 2')
plt.legend()
plt.savefig('t_sne_clusters.png')
plt.show()