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

def get_number_of_features(feature_vector):
    return sum(1 for x in feature_vector if x not in [-1, -1.0]) // 2


class Feature_cluster:
    def __init__(self, feature_vector, data_vector, t_sne):
       self.feature_vector = feature_vector
       self.data_vector = data_vector
       self.t_sne = t_sne 


file_list = ['parameters_all_a=0.8.json', 'parameters_all_a=0.9.json','parameters_all_a=0.95.json','parameters_all_a=1.0.json','parameters_all_a=1.05.json','parameters_all_a=1.1.json','parameters_all_a=1.2.json','parameters_all_a=1.5.json']


feature_vectors = []
lc_vectors = []
data = []

for file_path in file_list:
  with open(file_path, 'r') as f:
      datum = json.load(f)

  data.extend(datum)      
  # Extract feature vectors from the JSON structure
  feature_vectors_per_file = [item['feature_vector'] for item in datum]
  lc_vector_per_file = [item['light_curve'] for item in datum]
  
  print("opened "+ file_path + " with " + str(len(datum)) + " events")
  lc_vectors.extend(lc_vector_per_file)
  feature_vectors.extend(feature_vectors_per_file)



# Convert to DataFrame for easier handling
df = pd.DataFrame(feature_vectors)

# Fill NaN values with the column mean
#df.fillna(df.mean(), inplace=True)

lc_df = pd.DataFrame(lc_vectors)

# Extract feature matrix
feature_matrix = df.values
lc_matrix = lc_df.values


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

## Plot the t-SNE results with clusters
#plt.figure(figsize=(10, 6))
#colors = plt.cm.get_cmap('tab10', num_clusters)
#for cluster_num in range(num_clusters):
#    indices = np.where(clusters == cluster_num)[0]
#    plt.scatter(tsne_results[indices, 0], tsne_results[indices, 1], 
#                c=[colors(cluster_num)], s=50, alpha=0.7, label=f'Cluster {cluster_num}')
    

    
#plt.title('t-SNE Visualization of Feature Vectors with Clustering')
#plt.xlabel('t-SNE Component 1')
#plt.ylabel('t-SNE Component 2')
#plt.legend()
#plt.savefig('t_sne_clusters.png')
#plt.show()


## Plot the t-SNE results with num-of-features
#plt.figure(figsize=(10, 6))
#color_values = [get_number_of_features(feature_vectors[idx]) for idx in range(len(feature_vectors))]
#
#print(color_values)
#
##colors = plt.cm.get_cmap('tab10', len(color_values))
#plt.scatter(tsne_results[:, 0], tsne_results[:, 1], 
#                c=color_values, cmap='viridis',  s=50, alpha=0.7, label='color-coded')
#    
#plt.colorbar()
#
#plt.title('t-SNE Visualization of Feature Vectors with Clustering')
#plt.xlabel('t-SNE Component 1')
#plt.ylabel('t-SNE Component 2')
#plt.legend()
#plt.savefig('t_sne_num_of_features.png')
#plt.show()


feature_clusters = {}
data_clusters = {}
cluster_dict = {}
num_of_subclusters_dict = {}

for idx in range(len(feature_vectors)):
  num_of_features = get_number_of_features(feature_vectors[idx])  
  if num_of_features not in feature_clusters:
    feature_clusters[num_of_features] = [feature_vectors[idx]]
    data_clusters[num_of_features] = [data[idx]]
  else:
    feature_clusters[num_of_features].append(feature_vectors[idx])
    data_clusters[num_of_features].append(data[idx])

num_of_subcluster_dict = {1: 3, 2: 6, 3: 7, 4: 1, 5: 1}



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
  plt.savefig('t_sne_clusters_'+str(num_of_features)+'.png')  

  # Output the cluster results to separate files
  for cluster_num in range(num_clusters):
    cluster_indices = np.where(clusters == cluster_num)[0]
    cluster_data = [data_clusters[num_of_features][idx] for idx in cluster_indices]
    with open(os.path.join(output_dir, f'subcluster_{num_of_features}_{cluster_num}.json'), 'w') as f:
      json.dump(cluster_data, f, indent=2)

    cluster_params = [[item['a'], item['b'], item['theta'], item['q'], item['alpha']] for item in cluster_data]
    parameter_cluster_df = pd.DataFrame(cluster_params)
    if len(cluster_params) > 30:
      perplexity_params = min(30, len(cluster_params)-1)
      tsne_params = TSNE(n_components=2, random_state=42, perplexity = perplexity_params)
      tsne_params_results = tsne.fit_transform(parameter_cluster_df.values)
      plt.figure(figsize=(10, 6))      
      plt.scatter(tsne_params_results[:, 0], tsne_params_results[:, 1], 
                    c=[colors(cluster_num)], s=50, alpha=0.7, label=f'Cluster {num_of_features}_{cluster_num}')

      plt.title('t-SNE of parameters of lc subclusters')
      plt.xlabel('t-SNE Component 1')
      plt.ylabel('t-SNE Component 2')
      plt.legend()
      plt.savefig('t_sne_parameter_subclusters_'+str(num_of_features)+'_'+str(cluster_num)+'.png')  
    else:
        print('too little samples in cluster'+str(num_of_features)+'_'+str(cluster_num)+' : '+ str(len(cluster_params)))


  
      

     