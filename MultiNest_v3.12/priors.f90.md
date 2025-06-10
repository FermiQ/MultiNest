# Documentation for `priors.f90`

## Overview

The `priors.f90` module provides a collection of functions for transforming a uniformly distributed random variable (typically in the range [0,1]) into a sample from a specified prior distribution. These functions are essential for the nested sampling algorithm, as MultiNest internally works with points in a unit hypercube. The user's log-likelihood function is responsible for using these (or custom) transformations to map these unit hypercube points to physical parameter values before calculating the likelihood.

## Key Components

Each function takes a uniform random deviate `r` (double precision, nominally in [0,1]) and one or two parameters (`x1`, `x2`, `mu`, `sigma`, etc.) defining the shape of the target prior distribution.

### `DeltaFunctionPrior`
*   **Type**: Function
*   **Description**: Transforms `r` to a delta function at `x1`. Effectively, it always returns `x1`.
*   **Parameters**:
    *   `r`: (Input, `double precision`) Uniform random deviate.
    *   `x1`: (Input, `double precision`) The position of the delta function.
    *   `x2`: (Input, `double precision`) Unused.
*   **Returns**: `x1`.

### `UniformPrior`
*   **Type**: Function
*   **Description**: Transforms `r` to a sample from a uniform distribution between `x1` and `x2`.
*   **Parameters**:
    *   `r`: (Input, `double precision`) Uniform random deviate.
    *   `x1`: (Input, `double precision`) Lower bound of the uniform distribution.
    *   `x2`: (Input, `double precision`) Upper bound of the uniform distribution.
*   **Returns**: `x1 + r * (x2 - x1)`.

### `LogPrior`
*   **Type**: Function
*   **Description**: Transforms `r` to a sample from a log-uniform distribution between `x1` and `x2`. This means that `log10(parameter)` is uniformly distributed.
*   **Parameters**:
    *   `r`: (Input, `double precision`) Uniform random deviate.
    *   `x1`: (Input, `double precision`) Lower bound (positive).
    *   `x2`: (Input, `double precision`) Upper bound (positive).
*   **Returns**: `10**(log10(x1) + r * (log10(x2) - log10(x1)))`.

### `SinPrior`
*   **Type**: Function
*   **Description**: Transforms `r` to a sample from a distribution proportional to `sin(angle)`, where `angle` is in degrees between `x1` and `x2`. This is often used for inclination angles. The output is in radians.
*   **Parameters**:
    *   `r`: (Input, `double precision`) Uniform random deviate.
    *   `x1`: (Input, `double precision`) Lower bound of angle in degrees.
    *   `x2`: (Input, `double precision`) Upper bound of angle in degrees.
*   **Returns**: `acos(cos(x1*deg2rad) + r * (cos(x2*deg2rad) - cos(x1*deg2rad)))`.

### `CauchyPrior`
*   **Type**: Function
*   **Description**: Transforms `r` to a sample from a Cauchy (Lorentzian) distribution.
*   **Parameters**:
    *   `r`: (Input, `double precision`) Uniform random deviate.
    *   `x0`: (Input, `double precision`) Center (mean/median/mode) of the Cauchy distribution.
    *   `gamma`: (Input, `double precision`) Half-width at half-maximum (HWHM).
*   **Returns**: `x0 + gamma * tan(PI * (r - 0.5))`.

### `GaussianPrior`
*   **Type**: Function
*   **Description**: Transforms `r` to a sample from a Gaussian (normal) distribution using the inverse error function (`dierfc`).
*   **Parameters**:
    *   `r`: (Input, `double precision`) Uniform random deviate.
    *   `mu`: (Input, `double precision`) Mean of the Gaussian distribution.
    *   `sigma`: (Input, `double precision`) Standard deviation of the Gaussian distribution.
*   **Returns**: `mu + sigma * SqrtTwo * dierfc(2.0 * (1.0 - r))`.

### `LogNormalPrior`
*   **Type**: Function
*   **Description**: Transforms `r` to a sample from a log-normal distribution.
*   **Parameters**:
    *   `r`: (Input, `double precision`) Uniform random deviate.
    *   `a`: (Input, `double precision`) Mode of the log-normal distribution.
    *   `sigma`: (Input, `double precision`) Width parameter.
*   **Returns**: `a * dexp(sigma*sigma + sigma*SqrtTwo*dierfc(2.0*r))`.

### `ExponentialPrior`
*   **Type**: Function
*   **Description**: Transforms `r` to a sample from an exponential distribution.
*   **Parameters**:
    *   `r`: (Input, `double precision`) Uniform random deviate.
    *   `lambda`: (Input, `double precision`) Rate parameter.
*   **Returns**: `-log(r) / lambda`.

### `dierfc`
*   **Type**: Function
*   **Description**: Calculates the inverse of the complementary error function in double precision. This is a utility function used by `GaussianPrior` and `LogNormalPrior`.
*   **Parameters**:
    *   `y`: (Input, `double precision`) Value for which the inverse complementary error function is computed.
*   **Returns**: The inverse complementary error function of `y`.

## Important Variables/Constants

This module does not define global constants critical for external users. The behavior is primarily controlled by the parameters passed to its functions.

## Usage Examples

These functions are intended to be used within the user-supplied log-likelihood function to transform the unit hypercube parameters (elements of the `Cube` array, which are between 0 and 1) into physically meaningful parameter values.

```fortran
subroutine my_loglikelihood(Cube, n_dim, nPar, lnew, context_pass)
  use priors
  implicit none
  integer, intent(in) :: n_dim, nPar, context_pass
  double precision, intent(inout) :: Cube(nPar) ! Input: [0,1], Output: physical values
  double precision, intent(out) :: lnew

  double precision :: param1, param2, param3

  ! Example transformations:
  ! param1: Uniform prior between 10 and 20
  param1 = UniformPrior(Cube(1), 10.0d0, 20.0d0)
  Cube(1) = param1 ! Store physical value back if needed by dumper

  ! param2: Log-uniform prior between 1e-3 and 1e3
  param2 = LogPrior(Cube(2), 1.0d-3, 1.0d3)
  Cube(2) = param2

  ! param3: Gaussian prior with mean 0 and sigma 1
  param3 = GaussianPrior(Cube(3), 0.0d0, 1.0d0)
  Cube(3) = param3

  ! ... use param1, param2, param3 to calculate lnew ...
  lnew = -0.5d0 * (param1**2 + param2**2 + param3**2) ! Example
end subroutine my_loglikelihood
```

## Dependencies and Interactions

*   **`utils1.f90`**: The `dierfc` function is defined within this module but used by `GaussianPrior` and `LogNormalPrior`. (Correction: `dierfc` is actually in `priors.f90` itself).
*   **User's Log-Likelihood Function**: This is where these prior transformation functions are primarily used.
*   **`nested.F90`**: The `nested.F90` module provides the `Cube` array (with values in [0,1]) to the log-likelihood function.
```
