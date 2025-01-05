import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import os
import json
from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import SpectralClustering
#import hdbscan

from examples.analytic_classes import LightCurveAnalyzer

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

def get_number_of_features(feature_vector):
    return sum(1 for x in feature_vector[3:] if x not in [-1, -1.0]) // 3

def get_number_of_peaks(feature_vector):
    feature_vector_cut = feature_vector[3:]
    sum_of_peaks = 0
    for i in range(0, len(feature_vector_cut) // 3):
      if feature_vector_cut[3*i+2] > 0:
        sum_of_peaks += feature_vector_cut[3*i+2]
      else:
        break
    return sum_of_peaks


class Feature_cluster:
    def __init__(self, feature_vector, data_vector, t_sne):
       self.feature_vector = feature_vector
       self.data_vector = data_vector
       self.t_sne = t_sne 

file_path = "./ml_project/parameters_backup/32/"
file_list = ['parameters_a=0.75.json',
             'parameters_a=0.8.json',
             'parameters_a=0.85.json',
             'parameters_a=0.9.json',
             'parameters_a=0.95.json',
             'parameters_a=0.97.json',
             'parameters_a=0.98.json',
             'parameters_a=0.99.json',
             'parameters_a=1.0.json',
             'parameters_a=1.01.json',
             'parameters_a=1.02.json',
             'parameters_a=1.03.json',
             'parameters_a=1.05.json',
             'parameters_a=1.07.json',
             'parameters_a=1.1.json',
             'parameters_a=1.11.json',
             'parameters_a=1.12.json',
             'parameters_a=1.15.json',
             'parameters_a=1.2.json',
             'parameters_a=1.25.json',
             'parameters_a=1.3.json',
             'parameters_a=1.35.json',
             'parameters_a=1.4.json',
             'parameters_a=1.5.json']


feature_vectors = []
lc_vectors = []
data = []

for file in file_list:
  with open(file_path+file, 'r') as f:
      datum = json.load(f)

  # t-sne cannot work with NaNs, it is better to remove data that have NaNs in them directly
  clean_datum = [
      item for item in datum 
      if not (np.any(np.isnan(item['feature_vector'])) or np.any(np.isnan(item['light_curve'])) or get_number_of_features(item['feature_vector']) == 0)
  ]

  data.extend(clean_datum)      
  # Extract feature vectors from the JSON structure
  feature_vectors_per_file = [item['feature_vector'][:28] for item in clean_datum]
  lc_vector_per_file = [item['light_curve'] for item in clean_datum]
  
  print("opened "+ file_path + file + " with " + str(len(clean_datum)) + " events")
  lc_vectors.extend(lc_vector_per_file)
  feature_vectors.extend(feature_vectors_per_file)

# Convert to DataFrame for easier handling
df = pd.DataFrame(feature_vectors)

# Fill NaN values with the column mean
#df.fillna(-1.0, inplace=True)

lc_df = pd.DataFrame(lc_vectors)

# Extract feature matrix
feature_matrix = df.values
lc_matrix = lc_df.values


# Perform t-SNE
tsne = TSNE(n_components=2, random_state=42)
tsne_results = tsne.fit_transform(feature_matrix)

num_clusters = 7  # You can change this number based on your needs

agglomerative = AgglomerativeClustering(n_clusters=num_clusters, linkage='single')
clusters = agglomerative.fit_predict(tsne_results)

# Create a directory for cluster output files
cluster_output_dir = os.path.join('ml_project', 'cluster_results')
os.makedirs(cluster_output_dir, exist_ok=True)

fig_output_dir = './ml_project/'

# Plot the t-SNE results with clusters
plt.figure(figsize=(10, 6))
colors = plt.colormaps['tab10']
for cluster_num in range(len(clusters)):
    indices = np.where(clusters == cluster_num)[0]
    plt.scatter(tsne_results[indices, 0], tsne_results[indices, 1], 
                c=[colors(cluster_num)], s=50, alpha=0.7, label=f'Cluster {cluster_num}')
   
   
plt.title('t-SNE Visualization of Feature Vectors with Clustering')
plt.xlabel('t-SNE Component 1')
plt.ylabel('t-SNE Component 2')
plt.legend()
plt.savefig(fig_output_dir+'t_sne_clusters_all.png')
plt.close()


# Plot the t-SNE results with num-of-features
plt.figure(figsize=(10, 6))
color_values = [get_number_of_features(feature_vectors[idx]) for idx in range(len(feature_vectors))]
plt.scatter(tsne_results[:, 0], tsne_results[:, 1], 
                c=color_values, cmap='viridis',  s=50, alpha=0.7, label='color-coded')
plt.colorbar()
plt.title('t-SNE Visualization of Feature Vectors with Clustering')
plt.xlabel('t-SNE Component 1')
plt.ylabel('t-SNE Component 2')
plt.legend()
plt.savefig(fig_output_dir+'t_sne_num_of_features.png')
plt.close()

# Plot the t-SNE results with num-of-features
plt.figure(figsize=(10, 6))
color_values = [get_number_of_peaks(feature_vectors[idx]) for idx in range(len(feature_vectors))]
plt.scatter(tsne_results[:, 0], tsne_results[:, 1], 
                c=color_values, cmap='viridis',  s=50, alpha=0.7, label='color-coded')
plt.colorbar()
plt.title('t-SNE Visualization of Feature Vectors with Clustering')
plt.xlabel('t-SNE Component 1')
plt.ylabel('t-SNE Component 2')
plt.legend()
plt.savefig(fig_output_dir+'t_sne_num_of_peaks.png')
plt.close()


feature_clusters = {}
data_clusters = {}
cluster_dict = {}
num_of_subclusters_dict = {}

# Fill in dictionaries num_of_features : data
for idx in range(len(feature_vectors)):
  num_of_features = get_number_of_features(feature_vectors[idx])  
  if num_of_features not in feature_clusters:
    feature_clusters[num_of_features] = [feature_vectors[idx]]
    data_clusters[num_of_features] = [data[idx]]
  else:
    feature_clusters[num_of_features].append(feature_vectors[idx])
    data_clusters[num_of_features].append(data[idx])

num_of_subcluster_dict = {1: 6, 2: 4, 3: 11, 4: 6, 5: 1, 6: 1}

# Par cluster dict is a dictionary to represent number of sub-subclusters that emerge
# by clustering subcluster based on parameters
num_of_par_cluster_dict = {
  1: {  # Cluster 0
    0: 4,
    1: 4,
    2: 1,
    3: 1,
    4: 1,
    5: 1
  },
  2: {  # Cluster 0
    0: 2,
    1: 4,
    2: 1,
    3: 2
  },
  3: {  # Cluster 0
    0: 2,
    1: 1,
    2: 2,
    3: 2,
    4: 2,
    5: 1,
    6: 2,
    7: 1,
    8: 2,
    9: 1,
    10: 1
  },
  4: {  # Cluster 0
    0: 1,
    1: 1,
    2: 1,
    3: 1
  }, 
  5: {  # Cluster 0
    0: 2
  },
  6: {
    0: 1
  }    
}

for num_of_features, feature_cluster in feature_clusters.items():
  print("num_of_features: " +str(num_of_features) + ", num of events: " + str(len(feature_cluster))) 
  # Perform t-SNE
  feature_cluster_df = pd.DataFrame(feature_cluster)
  perplexity = min(30, len(feature_cluster)-1)
  tsne = TSNE(n_components=2, random_state=42, perplexity = perplexity)
  tsne_results = tsne.fit_transform(feature_cluster_df.values)

  if num_of_features in num_of_subcluster_dict:
    num_clusters = num_of_subcluster_dict[num_of_features]
  else:
    num_clusters = 1  # You can change this number based on your needs

  agglomerative = AgglomerativeClustering(n_clusters=num_clusters, linkage='single')
  clusters = agglomerative.fit_predict(tsne_results)

  plt.figure(figsize=(10, 6))
  colors = plt.colormaps['tab10']
  
  for cluster_num in range(num_clusters):
      indices = np.where(clusters == cluster_num)[0]
      plt.scatter(tsne_results[indices, 0], tsne_results[indices, 1], 
                  c=[colors(cluster_num)], s=50, alpha=0.7, label=f'Cluster {cluster_num}')

  plt.title('t-SNE Visualization of Feature Vectors with Clustering')
  plt.xlabel('t-SNE Component 1')
  plt.ylabel('t-SNE Component 2')
  plt.legend()
  plt.savefig(fig_output_dir+'t_sne_clusters_'+str(num_of_features)+'.png')  



  # For each subcluster, try to cluster even more based on the parameters
  for cluster_num in range(num_clusters):
    cluster_indices = np.where(clusters == cluster_num)[0]
    cluster_data = [data_clusters[num_of_features][idx] for idx in cluster_indices]
    with open(os.path.join(cluster_output_dir, f'subcluster_{num_of_features}_{cluster_num}.json'), 'w') as f:
      json.dump(cluster_data, f, indent=2)

    cluster_params = [[item['a'], item['b'], item['theta'],item['m2'], item['q'], item['alpha']] for item in cluster_data]
    parameter_cluster_df = pd.DataFrame(cluster_params)

    num_parclusters = 1

    if num_of_features in num_of_subcluster_dict:
      if cluster_num in num_of_par_cluster_dict[num_of_features]:
        num_parclusters = num_of_par_cluster_dict[num_of_features][cluster_num]

    if len(cluster_params) > 30:
      perplexity_params = min(30, len(cluster_params)-1)
      tsne_params = TSNE(n_components=2, random_state=42, perplexity = perplexity_params)
      tsne_params_results = tsne.fit_transform(parameter_cluster_df.values)

      paragglomerative = AgglomerativeClustering(n_clusters=num_parclusters, linkage='single')
      parclusters = paragglomerative.fit_predict(tsne_params_results)

      plt.figure(figsize=(10, 6))
      colors = plt.colormaps['tab10']

      plt.figure(figsize=(10, 6))

      for parcluster_num in range(num_parclusters):
        indices = np.where(parclusters == parcluster_num)[0]
        plt.scatter(tsne_params_results[indices, 0], tsne_params_results[indices, 1], 
                      c=[colors(parcluster_num)], s=50, alpha=0.7, label=f'Cluster {parcluster_num}')

      plt.title('t-SNE of parameters of lc subclusters ('+str(num_of_features)+'_'+str(cluster_num)+')')
      plt.xlabel('t-SNE Component 1')
      plt.ylabel('t-SNE Component 2')
      plt.legend()
      plt.savefig(fig_output_dir+'t_sne_parameter_subclusters_'+str(num_of_features)+'_'+str(cluster_num)+'.png')
      plt.close()  
    else:
        print('too little samples in cluster'+str(num_of_features)+'_'+str(cluster_num)+' : '+ str(len(cluster_params)))

    for parcluster_num in range(num_parclusters):
      parcluster_indices = np.where(clusters == parcluster_num)[0]
      parcluster_data = [data_clusters[num_of_features][idx] for idx in parcluster_indices]
      with open(os.path.join(cluster_output_dir, f'parcluster_{num_of_features}_{cluster_num}_{parcluster_num}.json'), 'w') as f:
        json.dump(parcluster_data, f, indent=2)
  
      

     