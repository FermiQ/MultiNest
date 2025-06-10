# Documentation for `utils1.f90`

## Overview

The `utils1.f90` module is a utility library for the MultiNest algorithm. It provides a wide range of helper functions primarily focused on numerical computations, matrix algebra, ellipsoid manipulation, and statistical calculations. These routines are used extensively throughout the other MultiNest modules.

## Key Components

### Matrix and Vector Operations
*   **`Diagonalize`**:
    *   Diagonalizes a symmetric matrix `a`. Returns eigenvalues in `diag` and eigenvectors in `a` (overwrites input). Uses LAPACK's `DSYEVR`.
    *   Manages workspace allocation for `DSYEVR` efficiently.
*   **`calc_covmat`**:
    *   Calculates the covariance matrix for a given set of points `p` and their `mean`.
*   **`calc_covmat_wt`**:
    *   Calculates a weighted covariance matrix.
*   **`calc_invcovmat`**:
    *   Calculates the inverse of a covariance matrix given its `evec` (eigenvectors) and `eval` (eigenvalues) as `inv_cov = evec * inv(diag(eval)) * evec^T`.
*   **`getTMatrix`**:
    *   Computes the transformation matrix `T` from eigenvectors and eigenvalues, such that `T * T^T = CovarianceMatrix`. Specifically, `TMatrix(i,j) = evec(i,j) * sqrt(eval(j))`, then transposed.

### Ellipsoid Operations and Point Properties
*   **`ScaleFactor`**:
    *   Calculates the maximum Mahalanobis radius (`kmax`) required to enclose all points `pt` with respect to a Gaussian distribution defined by `mean` and `inv_cov`. `kmax = max_i((pt_i - mean)^T * inv_cov * (pt_i - mean))`.
*   **`ptScaleFac`**:
    *   Calculates the Mahalanobis radius for a single point: `(pt - meanx)^T * invcovx * (pt - meanx)`.
*   **`MahaDis`**:
    *   Calculates the Mahalanobis distance of a point `pt` from an ellipsoid defined by `mean`, `inv_cov`, and an overall enlargement factor `kfac`. `MahaDis = ptScaleFac(...) / kfac`.
*   **`ellVol`**:
    *   Calculates the volume of an n-dimensional ellipsoid given its eigenvalues (`eval`) and an overall enlargement factor `k_fac`.
*   **`genPtInSpheroid`**:
    *   Generates a point uniformly distributed within an n-dimensional unit spheroid.
*   **`genPtOnSpheroid`**:
    *   Generates a point uniformly distributed on the surface of an n-dimensional unit spheroid.
*   **`genPtInEll`**:
    *   Generates a point uniformly distributed within a given ellipsoid (defined by `mean`, `efac` - enlargement factor, and `TMat` - transformation matrix).
*   **`genPtOnEll`**:
    *   Generates a point uniformly distributed on the surface of a given ellipsoid.
*   **`CalcEllProp`**:
    *   Calculates comprehensive properties of an ellipsoid bounding a set of points `pt`. Outputs include `mean`, `covmat`, `invcov`, `tMat`, `evec`, `eval`, `detcov`, point enlargement factor `kfac`, volume enlargement factor `eff`, and `vol`. It considers a minimum prior volume `pVol`.
*   **`CalcBEllInfo`**:
    *   Calculates properties for a bounding ellipsoid, similar to `CalcEllProp` but potentially with outlier removal and tailored for constructing bounding ellipsoids for clustering.
*   **`evolveEll`**:
    *   Updates ellipsoid properties (`kfac`, `eff`, `vol`) when a point is either added to or removed from the set defining the ellipsoid. It considers a target volume `pVol`.
*   **`enlargeEll`**:
    *   Enlarges an ellipsoid to encompass a `newpt`. Updates `kfac`, `eff`, `vol`, considering `pVol`.
*   **`ptIn1Ell`**:
    *   Logical function that checks if a given `pt` is inside an ellipsoid defined by `meanx`, `invcovx`, and `kfacx`.

### Mathematical and Statistical Utilities
*   **`LogSumExp`**:
    *   Calculates `log(exp(x) + exp(y))` robustly, avoiding underflow/overflow.
*   **`gammln`**: Computes `ln(Gamma(xx))`.
*   **`gammp`**: Computes the incomplete gamma function P(a,x).
*   **`gammq`**: Computes the complementary incomplete gamma function Q(a,x) = 1 - P(a,x).
*   **`gser`**: Series representation for `gammp`.
*   **`gcf`**: Continued fraction representation for `gammq`.
*   **`erf`**: Error function erf(x), computed using `gammp`.
*   **`stNormalCDF`**: Standard Normal Cumulative Distribution Function (CDF).
*   **`mNormalF`**: Calculates the multivariate normal probability density function value for a point `x` given `mean` and covariance matrix `C`. Uses Cholesky decomposition (`DPOTRF`).

### Miscellaneous Utilities
*   **`inprior`**:
    *   Logical function that checks if a point `p` (in hypercube coordinates) is within the unit hypercube [0,1]^np.
*   **`binSearch`**:
    *   Performs a binary search to find the position of `x` in a sorted integer `array`.
*   **`piksrt`**:
    *   Sorts an array `arr` and correspondingly rearranges a companion array `arr1`. (Insertion sort variant).
*   **`ceffEvolve`**:
    *   Evolves ellipsoid volume and enlargement factor in constant efficiency mode.

## Important Variables/Constants (Module Level)
*   `lwork`, `liwork`: (Integer) Store optimal workspace sizes for LAPACK routines, determined at runtime.
*   `setBlk`: (Logical) Flag to indicate if `lwork`/`liwork` have been set.

## Usage Examples

The routines in `utils1.f90` are fundamental building blocks for the core logic in `nested.F90`, `kmeans_clstr.f90`, `xmeans_clstr.f90`, and `posterior.F90`. They are not typically called directly by the end-user.

*   `CalcEllProp` and related functions are crucial for defining and sampling from ellipsoidal bounds in multimodal sampling.
*   `Diagonalize` and `calc_covmat` are used for Principal Component Analysis and characterizing point distributions.
*   `LogSumExp` is vital for numerically stable evidence calculations.

## Dependencies and Interactions

*   **`RandomNS` (from `utils.f90`)**: Used by `genPtInSpheroid` and `genPtOnSpheroid` for generating random numbers.
*   **LAPACK**: This module relies on LAPACK routines for linear algebra operations:
    *   `DSYEVR`: Used in `Diagonalize` for eigenvalue and eigenvector computation of a symmetric matrix.
    *   `DPOTRF`: Used in `mNormalF` for Cholesky decomposition.
*   **Other MultiNest Modules**: This module's functions are extensively called by:
    *   `nested.F90`: For ellipsoid calculations, random point generation within ellipsoids, matrix operations.
    *   `kmeans_clstr.f90`, `xmeans_clstr.f90`: For calculating cluster properties, distances, and diagonalizing covariance matrices.
    *   `posterior.F90`: For calculating multivariate normal probabilities or other statistical properties.
```
