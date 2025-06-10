# Documentation for `utils.f90`

## Overview

The `utils.f90` file defines a module named `RandomNS`. This module encapsulates a random number generator specifically tailored for use within the MultiNest algorithm. It provides routines for initializing the generator (potentially with multiple streams for parallel execution) and for generating both uniform and Gaussian random deviates.

The random number generator is based on a method proposed by George Marsaglia and was slightly modified by F. James. It's designed to produce long sequences of pseudo-random numbers, suitable for extensive Monte Carlo simulations.

## Key Components

### Module `RandomNS`

#### `initRandomNS`
*   **Type**: Subroutine
*   **Description**: Initializes the random number generator. It can set up multiple independent random number streams, which is useful for parallel (MPI) execution, ensuring that different processes use different sequences. If a seed `i` is provided, it's used to initialize the sequence; otherwise, system clock and date/time are used to generate initial seeds.
*   **Parameters**:
    *   `n`: (Input, `integer`) The number of random number streams to initialize (typically corresponding to the number of MPI nodes/threads).
    *   `i`: (Input, Optional, `integer`) An optional seed. If provided, it's used in conjunction with the stream index `k` to seed each stream. If not provided, seeds are derived from system time.

#### `killRandomNS`
*   **Type**: Subroutine
*   **Description**: Deallocates the memory used by the random number generator arrays (`C`, `CD`, `CM`, `U`, `I97`, `J97`, `ISET`, `GSET`). This should be called when the random number generator is no longer needed to free resources.

#### `RMARINNS`
*   **Type**: Subroutine
*   **Description**: This is an internal initialization routine for a specific stream of the random number generator. It sets up the internal state variables (`U` array, `C`, `CD`, `CM`, `I97`, `J97`) for the given stream `id` based on the input seeds `IJ` and `KL`.
*   **Parameters**:
    *   `IJ`: (Input, `integer`) The first seed value (0 <= IJ <= 31328).
    *   `KL`: (Input, `integer`) The second seed value (0 <= KL <= 30081).
    *   `id`: (Input, `integer`) The identifier for the random number stream to be initialized.

#### `ranmarns`
*   **Type**: Double Precision Function
*   **Description**: Generates a pseudo-random number uniformly distributed in the interval (0, 1). It uses the Marsaglia algorithm with its internal state arrays.
*   **Parameters**:
    *   `idg`: (Input, `integer`) The identifier for the random number stream to use (0-indexed, so `id = idg + 1` is used internally).
*   **Returns**: A `double precision` random number between 0 and 1.

#### `GAUSSIAN1NS`
*   **Type**: Double Precision Function
*   **Description**: Generates a pseudo-random number from a standard normal (Gaussian) distribution (mean 0, variance 1). It uses the Box-Muller transform method, utilizing two uniform random numbers from `ranmarns` to produce two Gaussian deviates (one is returned, the other is stored for the next call to improve efficiency).
*   **Parameters**:
    *   `idg`: (Input, `integer`) The identifier for the random number stream to use.
*   **Returns**: A `double precision` random number from a N(0,1) distribution.

## Important Variables/Constants (Module Level)

These are internal state variables for the random number generator and are not meant to be modified directly by the user.
*   `rand_instNS`: (Integer) Counter for instances, possibly to ensure different seeds if `initRandomNS` is called multiple times without an explicit seed.
*   `C`, `CD`, `CM`: (Double precision arrays) State variables for the Marsaglia generator.
*   `U`: (Double precision 2D array) Stores the state (recent random numbers) for each stream.
*   `I97`, `J97`: (Integer arrays) Pointers/indices for the `U` array for each stream.
*   `ISET`, `GSET`: (Integer and Double precision arrays) Used by `GAUSSIAN1NS` to store the second Gaussian deviate from the Box-Muller method.
*   `numNodes`: (Integer) Stores the number of initialized random number streams.

## Usage Examples

The `RandomNS` module is used internally by MultiNest whenever random numbers are required, for example:
*   Generating initial live points by sampling from the prior.
*   Generating trial points within ellipsoids during the sampling process.
*   Stochastic elements within clustering algorithms (like k-means initialization).

```fortran
! Conceptual internal usage
use RandomNS
implicit none

integer :: num_streams = 4, i
double precision :: u_rand, g_rand

! Initialize for 4 MPI processes (streams)
call initRandomNS(num_streams)

! On MPI rank 0 (idg = 0)
u_rand = ranmarns(0)
g_rand = GAUSSIAN1NS(0)

! ... later, cleanup ...
call killRandomNS()
```

## Dependencies and Interactions

*   This module is designed to be fairly standalone for random number generation.
*   It is used by various other MultiNest modules, including:
    *   `nested.F90` (e.g., in `getrandom`, `genPtInSpheroid` via `utils1`)
    *   `utils1.f90` (e.g., `genPtInSpheroid`, `genPtOnSpheroid`)
    *   `kmeans_clstr.f90` (for initializing k-means clusters)
```
