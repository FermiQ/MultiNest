# Documentation for `xmeans_clstr.f90`

## Overview

The `xmeans_clstr.f90` module implements the X-means clustering algorithm, a variation of K-means that automatically estimates the optimal number of clusters 'K'. It achieves this by repeatedly applying K-means (typically starting with K=2 for a parent cluster) and using a model selection criterion, often the Bayesian Information Criterion (BIC) or an Anderson-Darling test for normality, to decide whether to split a cluster into sub-clusters or keep the parent cluster.

In MultiNest, this module is likely used for more adaptive and automated mode identification during the nested sampling process, especially in multimodal scenarios.

## Key Components

### Module `xmeans_clstr`

The module provides several `doXmeans` and `doGmeans` wrapper subroutines that call the core recursive functions (`Xmeans`, `Gmeans`, `makeDino`).

#### `doXmeans`
*   **Type**: Subroutine
*   **Description**: A wrapper for the X-means algorithm. It initializes and calls the main `Xmeans` recursive routine. It rearranges points and their likelihoods according to the found clusters.
*   **Parameters**: `points`, `npt`, `np` (dimensionality), `nClstr` (output number of clusters), `ptClstr` (output points per cluster), `cMean`, `cCovmat`, etc. (output cluster properties), `like` (likelihoods), `min_pt`.

#### `doXmeans5`
*   **Type**: Subroutine
*   **Description**: X-means with shared points. Points are rearranged, and auxiliary variables (`auxa`) associated with points are also rearranged. It calls `Xmeans5`.
*   **Parameters**: `points`, `npt`, `np`, `nClstr`, `ptClstr`, `ad` (shared points), `naux` (num auxiliary vars), `auxa`, `min_pt`.

#### `doXmeans8`
*   **Type**: Subroutine
*   **Description**: X-means variant that takes an initial set of clusters and further subdivides them.
*   **Parameters**: `points`, `like`, `npt`, `np`, `nClstr` (input/output), `ptClstr` (input/output), `ad`, `min_pt`. Calls `Xmeans8`.

#### `doXmeans6` / `doXmeans7`
*   **Type**: Subroutine
*   **Description**: Simplified X-means wrappers, likely assuming isotropic Gaussians or using simpler splitting criteria. `doXmeans6` also handles an auxiliary array. They call `Xmeans6` and `Xmeans7` respectively.

#### `doGmeans` / `doGmeans2`
*   **Type**: Subroutine
*   **Description**: Wrappers for G-means clustering, which uses a statistical test (Anderson-Darling) to check if a cluster should be split. `doGmeans` handles an auxiliary array. They call `Gmeans` and `Gmeans2`.

#### `Xmeans` (Recursive Subroutine)
*   **Type**: Recursive Subroutine
*   **Description**: The core X-means logic. For a given set of points, it:
    1.  Calculates properties assuming one cluster (BIC1).
    2.  Splits the points into two sub-clusters using K-means (typically `kmeans3`).
    3.  Calculates properties assuming two clusters (BIC2).
    4.  If BIC1 > BIC2 (i.e., two clusters are preferred) and other conditions (min points per cluster, max total clusters) are met, it recursively calls `Xmeans` on the two sub-clusters.
    5.  Otherwise, it considers the current set of points as a single cluster and stores its properties.
*   **Parameters**: `pt`, `npt`, `like`, `min_pt`.

#### `Xmeans4` / `Xmeans5` / `Xmeans8` / `Xmeans6` / `Xmeans7` / `Xmeans9` (Recursive Subroutines)
*   **Type**: Recursive Subroutines
*   **Description**: Variants of the core X-means logic, differing in how they handle shared points (`Xmeans4`, `Xmeans5`, `Xmeans9`), auxiliary data (`Xmeans6`), initial cluster division (`Xmeans8`), or the splitting criterion (e.g., using simpler isotropic BIC in `Xmeans6`, `Xmeans7`).

#### `Gmeans` / `Gmeans2` (Recursive Subroutines)
*   **Type**: Recursive Subroutines
*   **Description**: Core G-means logic. Similar to `Xmeans`, but instead of BIC, it uses the Anderson-Darling test on the projection of data along the vector connecting the means of the two potential sub-clusters. If the projected data is not normally distributed, the split is accepted.
*   **Parameters**: `pt`, `npt`, `naux`, `auxa` (for `Gmeans`), `min_pt`.

#### `Dinosaur`
*   **Type**: Logical Function
*   **Description**: A more complex clustering/splitting routine. It seems to iteratively refine clusters using `makeDino` and `postDino`. It considers prior volumes, current volumes, and volume tolerances to decide on the final set of clusters.
*   **Parameters**: Extensive, including point data, auxiliary data, ellipsoid properties, prior/current volumes, and control parameters.
*   **Return Value**: `.true.` if a satisfactory clustering is found, `.false.` otherwise.

#### `makeDino` (Recursive Subroutine)
*   **Type**: Recursive Subroutine
*   **Description**: Part of the "Dinosaur" clustering. It attempts to split a given cluster into `nk` (typically 2) sub-clusters using `Dmeans` from `kmeans_clstr`. If the split is favorable (based on volume reduction or other criteria), it recursively calls itself on the sub-clusters.
*   **Parameters**: Includes current cluster's points, properties, and target prior volume for sub-clusters.

#### `postDino`
*   **Type**: Logical Function
*   **Description**: Post-processes clusters generated by `makeDino`. It may reassign points between clusters if a cluster's volume relative to its expected prior volume is too large (`fVal > fTol`). It can also merge clusters.
*   **Return Value**: `.true.` if further processing by `makeDino` might be needed (e.g., if a large, diffuse cluster remains).

#### `calcBIC1` / `calcBIC1_iso`
*   **Type**: Function
*   **Description**: Calculates the Bayesian Information Criterion (BIC) for a single cluster model. `calcBIC1_iso` assumes isotropic Gaussians.
*   **Parameters**: Point data, number of points, mean, inverse covariance, determinant of covariance (or variance for isotropic).

#### `calcBIC2` / `calcBIC2_iso`
*   **Type**: Function
*   **Description**: Calculates the BIC for a two-cluster model. `calcBIC2_iso` assumes isotropic Gaussians.
*   **Parameters**: Data for two clusters (points, counts, means, inv. covariances, det. covariances or variances).

#### `AndersonDarling`
*   **Type**: Logical Function
*   **Description**: Performs an Anderson-Darling test for normality on data points projected onto the vector connecting the means of two potential sub-clusters.
*   **Return Value**: `.true.` if the data is consistent with a single normal distribution (no split), `.false.` if the data suggests two distributions (split).

## Important Variables/Constants (Module Level)

*   `n_dim`: (Integer) Stores the dimensionality of the data being clustered.
*   `maxClstr`: (Integer) Maximum number of clusters allowed.
*   Various allocatable arrays (e.g., `xclsMean`, `xclsCovmat`, `ptInClstr`) to store properties of identified clusters during the recursive process.

## Usage Examples

The subroutines in `xmeans_clstr.f90` are primarily intended for internal use by the MultiNest algorithm, particularly within the `clusteredNest` routine in `nested.F90` when `multimodal` sampling is enabled. They provide a more automated way to find and separate modes compared to fixed-K k-means.

## Dependencies and Interactions

*   **`kmeans_clstr.f90`**:
    *   Uses `kmeans3` (simple k-means) as a subroutine for splitting clusters in `Xmeans` and `Gmeans`.
    *   The `Dmeans` function from `kmeans_clstr` is used within `makeDino`.
*   **`utils1.f90`**:
    *   Uses functions for matrix operations (`calc_covmat`, `diagonalize`, `calc_invcovmat`).
    *   Uses ellipsoid property calculation routines (`CalcEllProp`, `MahaDis`).
    *   Uses statistical tests/CDFs (`stNormalCDF`).
*   **`nested.F90`**: The X-means routines are called from the main nested sampling loop (`clusteredNest` via `Dinosaur`) to help identify and sample multiple modes in the parameter space.
```
