# Documentation for `kmeans_clstr.f90`

## Overview

This Fortran module, `kmeans_clstr`, provides a suite of k-means clustering algorithms. K-means is a method for partitioning a dataset into 'k' distinct, non-overlapping clusters. In the context of MultiNest, these algorithms are likely employed to identify modes (peaks) or distinct regions in the parameter space being explored by the nested sampling algorithm. This helps in efficiently sampling multimodal distributions.

The module includes several variations of k-means:
*   Soft k-means (axis-aligned and spherical Gaussians)
*   Simple k-means
*   Incremental K-means (including versions that can determine K)

## Key Components

### `kmeans`
*   **Type**: Subroutine
*   **Description**: Implements a soft k-means clustering algorithm. It models clusters as axis-aligned Gaussians.
*   **Parameters**:
    *   `k`: (Input) The desired number of clusters.
    *   `pt`: (Input) `double precision pt(numdim,npt)` - The data points to cluster.
    *   `npt`: (Input) Number of data points.
    *   `numdim`: (Input) Dimensionality of the data points.
    *   `cluster`: (Output) `integer cluster(npt)` - Array storing the cluster assignment for each point.
    *   `min_pt`: (Input) Minimum number of points allowed in a cluster.

### `kmeans2`
*   **Type**: Subroutine
*   **Description**: Implements a soft k-means clustering algorithm, modeling clusters as spherical Gaussians.
*   **Parameters**: Same as `kmeans`.

### `kmeans3`
*   **Type**: Subroutine
*   **Description**: Implements a simple, hard k-means clustering algorithm.
*   **Parameters**:
    *   `k`: (Input) The desired number of clusters.
    *   `pt`: (Input) `double precision pt(numdim,npt)` - The data points.
    *   `npt`: (Input) Number of data points.
    *   `numdim`: (Input) Dimensionality.
    *   `means`: (Output) `double precision means(k,numdim)` - Calculated means of the clusters.
    *   `cluster`: (Output) `integer cluster(npt)` - Cluster assignment for each point.
    *   `min_pt`: (Input) Minimum number of points allowed in a cluster.

### `kmeans4`
*   **Type**: Subroutine
*   **Description**: Implements an Incremental K-means algorithm based on Pham, Dimov, Nguyen (2005). This version starts with one cluster and incrementally adds new clusters.
*   **Parameters**: Same as `kmeans`.

### `kmeans5`
*   **Type**: Subroutine
*   **Description**: Implements an Incremental K-means algorithm with unknown K, also based on Pham, Dimov, Nguyen (2005). It uses an evaluation function to determine the optimal number of clusters.
*   **Parameters**:
    *   `k`: (Input/Output) Initial guess for K / Output determined K.
    *   `pt`, `npt`, `numdim`, `cluster`, `min_pt`: Same as `kmeans`.

### `kmeans6`
*   **Type**: Subroutine
*   **Description**: A 2-means algorithm (k=2) where clusters are initialized at their expected positions based on principal component analysis.
*   **Parameters**:
    *   `pt`, `npt`, `numdim`, `cluster`: Same as `kmeans`.

### `kmeans7`
*   **Type**: Subroutine
*   **Description**: An Incremental K-means algorithm. If the number of points in any cluster is less than `numdim+1` (or `min_p`), it borrows the closest points from its nearest neighbor cluster. It can also identify shared points between clusters.
*   **Parameters**:
    *   `k`: (Input/Output) Desired/Final number of clusters.
    *   `pt`: (Input) Data points.
    *   `npt`: (Input) Number of points.
    *   `numdim`: (Input) Dimensionality.
    *   `cluster`: (Output) Cluster assignments.
    *   `ad`: (Input) Number of shared points to find.
    *   `cluster2`: (Output) `integer cluster2(:,:)` - Information about shared points.
    *   `min_p`: (Input) Minimum points per cluster.

### `sKmeans` (Likely a simplified or specific version of `kmeans7`)
*   **Type**: Subroutine
*   **Description**: Appears very similar to `kmeans7`, implementing an incremental K-means that borrows points from neighboring clusters if a cluster size falls below `min_p`.
*   **Parameters**: Same as `kmeans7` (excluding `ad` and `cluster2`).

### `kmeans8`
*   **Type**: Subroutine
*   **Description**: Implements a soft k-means algorithm corresponding to a model of arbitrary Gaussians (non-axis-aligned).
*   **Parameters**:
    *   `k`, `pt`, `npt`, `numdim`, `cluster`, `min_pt`: Same as `kmeans`.
    *   `means`: (Output) `double precision means(k,numdim)` - Calculated means.
    *   `Covmat`: (Output) `double precision covmat(k,numdim,numdim)` - Calculated covariance matrices for each cluster.

### `Dmeans`
*   **Type**: Function
*   **Description**: "Dinosaur clustering." This function appears to be a more sophisticated clustering or mode splitting algorithm. It takes existing ellipsoid properties (means, covariances, etc.) of `k-1` clusters and attempts to find a better partition into `k` clusters, possibly by splitting an existing ellipsoid.
*   **Parameters**: Extensive list including properties of `k-1` clusters and outputs properties for `k` clusters. Key inputs include `pt` (points), `like` (scaled log-likelihoods), `ndim`, `min_pt`, `pVol` (prior volume).
*   **Return Value**: `integer Dmeans` - Status code (0 for better partition, 1 for partition found but not better, 2 if couldn't partition).

### `kmeansDis`
*   **Type**: Subroutine
*   **Description**: An Incremental K-means subroutine where a specific existing cluster `m` is targeted for splitting. It returns means and distortions.
*   **Parameters**:
    *   `k`: (Input) Current number of clusters.
    *   `pt`, `npt`, `numdim`, `cluster`: As in `kmeans`.
    *   `m`: (Input) The cluster index to be split.
    *   `totR`: (Input/Output) Total points in each cluster.
    *   `distortion`: (Output) Distortion of each cluster.
    *   `means`: (Input/Output) Cluster means.

## Important Variables/Constants

This module does not define global constants critical for external users. The behavior is primarily controlled by the parameters passed to its subroutines.

## Usage Examples

The subroutines in `kmeans_clstr.f90` are primarily intended for internal use by the MultiNest algorithm (e.g., within `nested.F90` or `xmeans_clstr.f90`) for tasks like identifying and separating modes in the parameter space. Direct usage by an end-user of MultiNest is unlikely.

Example (conceptual, how `kmeans3` might be called internally):
```fortran
! integer, parameter :: num_dimensions = 3, num_points = 100, k_clusters = 4, min_pts_in_cluster = 5
! double precision :: points_to_cluster(num_dimensions, num_points)
! double precision :: cluster_means(k_clusters, num_dimensions)
! integer :: point_assignments(num_points)

! ! ... (populate points_to_cluster) ...

! call kmeans3(k_clusters, points_to_cluster, num_points, num_dimensions, &
!              cluster_means, point_assignments, min_pts_in_cluster)
```

## Dependencies and Interactions

*   **`RandomNS` (from `utils.f90`)**: Used for random number generation, essential for initializing cluster means or other stochastic parts of the algorithms.
*   **`utils1.f90`**: Likely uses utility functions from this module for tasks such as:
    *   Calculating covariance matrices (`calc_covmat`).
    *   Diagonalizing matrices (`diagonalize`).
    *   Calculating Mahalanobis distances or multivariate normal probabilities (`mNormalF`).
*   **`xmeans_clstr.f90`**: Some k-means variants here might be used as components within the X-means algorithm.
*   **`nested.F90`**: The clustering algorithms are called from the main nested sampling routines to help with mode separation and efficient exploration of complex parameter spaces.
```
