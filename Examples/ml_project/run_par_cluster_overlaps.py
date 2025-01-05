# this script reads number of parameter files clustered according to lens even parameters
# and looks if the clusters have overlaps and then lists them.

#!/usr/bin/env python3

import os
import glob
import json
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

# As an alternative to alpha-shape 
import alphashape
import pyvista as pv

##############################################################################
# 1. LOAD CLUSTERS FROM JSON FILES (3D: a, q, alpha)
##############################################################################
def load_clusters_from_directory(directory_path):
    """
    Reads all JSON files named 'parcluster_i_j_k.json' in the given directory.
    Extracts the 3D parameter vector [a, q, alpha] for each point in the cluster.
    Returns a dictionary:

        {
            "i_j_k": np.array([
                [a1, q1, alpha1],
                [a2, q2, alpha2],
                ...
            ]),
            ...
        }

    where "i_j_k" is the cluster identifier string, e.g. "1_2_3".
    """
    cluster_dict = {}

    # Pattern to match parcluster_i_j_k.json
    search_pattern = os.path.join(directory_path, "parcluster_*_*_*.json")
    json_files = glob.glob(search_pattern)

    for filename in json_files:
        # Extract the cluster ID from the filename, e.g. parcluster_1_2_3.json -> "1_2_3"
        base_name = os.path.basename(filename)
        cluster_id = base_name.replace("parcluster_", "").replace(".json", "")

        # Parse JSON file (each file is an array of dicts)
        with open(filename, "r") as f:
            data_list = json.load(f)

        # Gather the relevant 3D data: a, q, alpha
        points = []
        for entry in data_list:
            a     = entry["a"]
            q     = entry["q"]
            alpha = entry["alpha"]
            points.append([a, q, alpha])

        cluster_dict[cluster_id] = np.array(points)

    return cluster_dict


##############################################################################
# 2. BUILD A CONVEX HULL IN 3D
##############################################################################
def build_convex_hull(points_3d):
    """
    Given an (N x 3) array 'points_3d', build a 3-dimensional ConvexHull
    using scipy.spatial.ConvexHull. Returns the ConvexHull object or None
    if degenerate.
    """
    if points_3d is None or len(points_3d) < 4:
        # In 3D, we typically need at least 4 non-coplanar points for a volume.
        # But we can attempt with fewer, though it might have zero volume.
        # Let's just try it. If Qhull fails, return None.
        pass

    try:
        hull = ConvexHull(points_3d)
        return hull
    except:
        # Qhull can fail for degenerate inputs (coplanar, collinear, etc.)
        return None


##############################################################################
# 3. MONTE CARLO APPROXIMATION OF INTERSECTION VOLUME (3D)
##############################################################################
def intersection_volume_monte_carlo(hullA, hullB, pointsA, pointsB, num_samples=5000):
    """
    Approximate the intersection volume of two convex hulls in 3D by random sampling:
      1) Build an axis-aligned bounding box from pointsA + pointsB.
      2) Draw 'num_samples' random points in that bounding box.
      3) Check how many are inside both hulls (via a naive 'point in hull' check).
      4) Multiply that fraction by the bounding box volume to estimate intersection.

    Returns an approximate intersection volume in 3D.
    """
    if hullA is None or hullB is None:
        return 0.0

    # Combine points to get bounding box
    all_points = np.vstack((pointsA, pointsB))
    mins = np.min(all_points, axis=0)
    maxs = np.max(all_points, axis=0)
    box_volume = np.prod(maxs - mins)

    if box_volume <= 0:
        return 0.0

    volA = hullA.volume  # 3D volume of hull A
    volB = hullB.volume  # 3D volume of hull B
    if volA == 0 or volB == 0:
        return 0.0

    # Extract the points used for hull-vertices
    hullA_points = pointsA[hullA.vertices]
    hullB_points = pointsB[hullB.vertices]

    rng = np.random.default_rng()
    random_points = rng.uniform(mins, maxs, size=(num_samples, 3))

    inside_count = 0
    for pt in random_points:
        if (is_point_in_hull(pt, hullA_points, volA) and
            is_point_in_hull(pt, hullB_points, volB)):
            inside_count += 1

    fraction = inside_count / num_samples
    return fraction * box_volume


def is_point_in_hull(pt, hull_points, hull_volume, tol=1e-12):
    """
    Naive 3D "point in hull" check via volume comparison:
      - Build hull with (hull_points + [pt])
      - If new hull's volume == old hull's volume (within tolerance), 'pt' is inside.

    This is repeated for each test point, so it can be slow. 
    For larger data, consider a half-space approach instead.
    """
    # If hull_points is too small, Qhull won't form a volume.
    if len(hull_points) < 4:
        return False

    test_set = np.vstack([hull_points, pt])
    try:
        test_hull = ConvexHull(test_set)
        return np.isclose(test_hull.volume, hull_volume, atol=tol)
    except:
        return False


##############################################################################
# 4. MAIN SCRIPT
##############################################################################
def main(directory="ml_project/cluster_results",
         output_csv="intersection_volumes_3d.csv",
         num_samples=5000):
    """
    1) Load cluster data from JSON files in 'directory' (keys: a, q, alpha).
    2) Build 3D hull for each cluster.
    3) Compute pairwise intersection volumes (via Monte Carlo).
    4) Output to CSV a table with rows/columns as cluster IDs, and
       each cell = intersection volume in 3D.
    """
    # Step 1: Load cluster data => dict { cluster_id : Nx3 array }
    cluster_dict = load_clusters_from_directory(directory)

    # Sort cluster IDs for consistent ordering in the output
    cluster_ids = sorted(cluster_dict.keys())
    n = len(cluster_ids)

    # Prepare a 2D matrix of intersection volumes
    intersection_matrix = np.zeros((n, n), dtype=float)

    ## Step 2: Build up dictionary of meshes 
    meshes = {}
    for cid, points_3d in cluster_dict.items():
        if points_3d.shape[0] < 4:
            print(f"'{cid}' not enough points to form tetrahedron'")
            meshes[cid] = None
            continue
        
        alpha = 1.0
        #mesh_alpha = alphashape.alphashape(points_3d, alpha)
        #mesh_wrap = pv.wrap(mesh_alpha)
        #mesh_wrap.clean()
        #meshes[cid] = pv.wrap(mesh_wrap)
        cloud = pv.PolyData(points_3d)
        tet_mesh = cloud.delaunay_3d(alpha)
        
        if tet_mesh.n_cells == 0:
            # It's empty -> skip or volume=0
            meshes[cid] = None
            print(f"'{cid}' tet has zero cells'")
            continue  
        
        surface_mesh = tet_mesh.extract_surface()
        # Force triangulation (make sure all faces are triangles)
        surface_mesh = surface_mesh.triangulate()

        # Optionally clean to remove any stray/duplicate cells
        surface_mesh = surface_mesh.clean()
        # remove 90% of sells
        surface_mesh = surface_mesh.decimate(0.8)

        if surface_mesh.n_cells == 0:
            # It's empty -> skip or volume=0
            meshes[cid] = None
            print(f"'{cid}' surface has zero cells'")
            continue 

        meshes[cid] = surface_mesh
        print(f"'{cid}' successful'")


    for i in range(n):
        cidA =  cluster_ids[i]
        for j in range(n):
            cidB = cluster_ids[j]
        
            if meshes[cidA] is None or meshes[cidB] is None:
                 intersection_matrix[i, j] = 0.0
                 print(f"'{cidA}' or '{cidB}' is None'")
                 continue

            if meshes[cidA].volume == 0.0 or meshes[cidB].volume == 0.0:
               intersection_matrix[i, j] = 0.0      
               print(f"'{cidA}' or '{cidB}' is zero volume'")
               continue

            if i == j:
                intersection_matrix[i, j] = meshes[cidA].volume
            else:
                try:
                    intersection = meshes[cidA].boolean_intersection(meshes[cidB], tolerance=1e-4)
                    if intersection is not None:
                        intersection_matrix[i, j] = intersection.volume
                    else:
                        print(f"Intersection of '{cidA}' and '{cidB}' None")
                except Exception as e:
                    # If the intersection fails for any reason, treat as zero overlap
                    print(f"Boolean intersection failed in Python with: {e}")
                    intersection_matrix[i, j] = 0.0

            print(f"Intersection of '{cidA}' and '{cidB}' is '{intersection_matrix[i, j]}'")

 
    # Step 3: Write the results to CSV
    df = pd.DataFrame(intersection_matrix, index=cluster_ids, columns=cluster_ids)
    df.to_csv(output_csv)
    print(f"Saved 3D intersection volume table to '{output_csv}'")


##############################################################################
# 5. ENTRY POINT
##############################################################################
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Compute pairwise intersection volumes in 3D parameter space (a,q,alpha)."
    )
    parser.add_argument("--directory", type=str, default="ml_project/cluster_results",
                        help="Path to directory containing parcluster_i_j_k.json files.")
    parser.add_argument("--output_csv", type=str, default="intersection_volumes_3d.csv",
                        help="Output CSV file for intersection volume table in 3D.")
    parser.add_argument("--num_samples", type=int, default=5000,
                        help="Number of Monte Carlo samples for intersection volume.")
    args = parser.parse_args()

    main(directory=args.directory,
         output_csv=args.output_csv,
         num_samples=args.num_samples)
