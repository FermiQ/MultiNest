# Documentation for `nested.F90`

## Overview

The `nested.F90` file is the heart of the MultiNest algorithm. It contains the primary module `Nested`, which orchestrates the entire nested sampling process. This includes initialization, management of live points, the main sampling loop, handling of MPI for parallel execution, and interfacing with user-provided log-likelihood and dumper functions (via the C wrapper if necessary).

## Key Components

### Module `Nested`

#### `nestRun`
*   **Type**: Subroutine
*   **Description**: This is the main public entry point for initiating a MultiNest run. It sets up global parameters, initializes MPI (if configured and `initMPI` is true), handles resume logic, and calls the internal `Nestsample` routine to perform the sampling.
*   **Parameters**:
    *   `nest_IS`: (Input, `logical`) Importance Sampling mode.
    *   `nest_mmodal`: (Input, `logical`) Multimodal sampling mode.
    *   `nest_ceff`: (Input, `logical`) Constant efficiency mode.
    *   `nest_nlive`: (Input, `integer`) Number of live points.
    *   `nest_tol`: (Input, `double precision`) Tolerance for evidence calculation. When `log(Z_i) - log(Z_{i-1}) < nest_tol`, the run terminates.
    *   `nest_ef`: (Input, `double precision`) Target efficiency for sampling.
    *   `nest_ndims`: (Input, `integer`) Number of dimensions (parameters being varied by MultiNest).
    *   `nest_totPar`: (Input, `integer`) Total number of parameters (including derived parameters passed to the log-likelihood).
    *   `nest_nCdims`: (Input, `integer`) Number of parameters to be used for clustering.
    *   `maxClst`: (Input, `integer`) Maximum number of modes for multimodal runs.
    *   `nest_updInt`: (Input, `integer`) Interval for updating output files.
    *   `nest_Ztol`: (Input, `double precision`) Tolerance for Z in posterior_sampling. Samples with `log(Z_local) < Ztol` are ignored.
    *   `nest_root`: (Input, `character(LEN=1000)`) Root for output file names.
    *   `seed`: (Input, `integer`) Random number seed. `<0` for system clock.
    *   `nest_pWrap`: (Input, `integer array`) Indicates if a parameter is periodic (wraparound).
    *   `nest_fb`: (Input, `logical`) Feedback to stdout on/off.
    *   `nest_resume`: (Input, `logical`) Resume from previous run if files exist.
    *   `nest_outfile`: (Input, `logical`) Write output files to disk.
    *   `initMPI`: (Input, `logical`) If true, initializes MPI. Set to false if MPI is already initialized externally.
    *   `nest_logZero`: (Input, `double precision`) Log-likelihood value assigned to points outside the prior.
    *   `nest_maxIter`: (Input, `integer`) Maximum number of iterations. `0` for no limit.
    *   `loglike`: (Input, Subroutine Interface) User-provided log-likelihood function.
    *   `dumper`: (Input, Subroutine Interface) User-provided dumper function (called at `nest_updInt` intervals).
    *   `context`: (Input, `integer`) Context variable passed to `loglike` and `dumper`.

#### `Nestsample`
*   **Type**: Subroutine
*   **Description**: This internal subroutine manages the core nested sampling logic after initial setup. It either generates initial live points (if not resuming or if live points weren't fully generated) or proceeds directly to the main sampling loop (`clusteredNest`).
*   **Parameters**:
    *   `loglike`, `dumper`, `context`: Same as in `nestRun`.

#### `gen_initial_live`
*   **Type**: Subroutine
*   **Description**: Generates the initial set of `nlive` live points by sampling from the prior distribution and calculating their likelihoods. It handles resuming live point generation if a previous run was interrupted. MPI is used to distribute the generation of initial live points.
*   **Parameters**:
    *   `p`: (Output, `double precision array`) Array to store hypercube coordinates of live points.
    *   `phyP`: (Output, `double precision array`) Array to store physical coordinates of live points.
    *   `l`: (Output, `double precision array`) Array to store log-likelihoods of live points.
    *   `loglike`, `dumper`, `context`: Same as in `nestRun`.

#### `getrandom`
*   **Type**: Subroutine
*   **Description**: A simple utility to generate a random point within the unit hypercube using `ranmarNS`.
*   **Parameters**:
    *   `n`: (Input, `integer`) Number of dimensions.
    *   `x`: (Output, `double precision array`) The generated random point in the unit hypercube.
    *   `id`: (Input, `integer`) MPI rank or thread ID for the random number generator.

#### `clusteredNest`
*   **Type**: Subroutine
*   **Description**: This is the main iterative loop of the nested sampling algorithm. In each iteration, it:
    1.  Identifies the live point with the lowest likelihood.
    2.  Adds this point to the set of dead points, contributing to the evidence calculation.
    3.  Generates a new live point by sampling from the prior within the likelihood constraint defined by the removed point. This step involves mode identification (clustering using k-means or X-means) and ellipsoidal sampling if `multimodal` or `ceff` modes are active.
    4.  Updates evidence and other statistics.
    5.  Checks for termination conditions.
*   **Parameters**:
    *   `p`, `phyP`, `l`, `l0`: (Input/Output) Arrays for live point hypercube coordinates, physical coordinates, log-likelihoods, and previous lowest log-likelihoods.
    *   `loglike`, `dumper`, `context`: Same as in `nestRun`.

## Important Variables/Constants (Module Level)

Many of these are set by `nestRun` based on its input parameters.
*   `my_rank`: (Integer) MPI rank of the current process.
*   `mpi_nthreads`: (Integer) Total number of MPI processors.
*   `nlive`: (Integer) Number of live points.
*   `ndims`: (Integer) Number of dimensions (parameters).
*   `totPar`: (Integer) Total number of parameters (physical).
*   `nCdims`: (Integer) Number of parameters for clustering.
*   `logZero`: (Double precision) Log-likelihood for points outside prior.
*   `maxIter`: (Integer) Maximum iterations.
*   `tol`: (Double precision) Evidence tolerance for termination.
*   `ef`: (Double precision) Target sampling efficiency.
*   `multimodal`: (Logical) True if multimodal sampling is enabled.
*   `ceff`: (Logical) True if constant efficiency mode is enabled.
*   `IS`: (Logical) True if Importance Sampling is enabled.
*   `outfile`: (Logical) True if output files should be written.
*   `resume`: (Logical) True if resuming from a previous run.
*   `broot`: (Character) Root for output file names.
*   `pWrap`: (Logical array) Indicates periodic parameters.

## Usage Examples

The primary way to use the `Nested` module is by calling the `nestRun` subroutine.

```fortran
! Example Fortran call (conceptual)
program my_analysis
  use Nested
  implicit none

  ! Define my_loglikelihood and my_dumper subroutines
  ! ...

  integer, parameter :: Ndims = 5, Nlive = 1000, TotPar = 5
  integer :: pWrap(Ndims)
  pWrap = 0 ! Assuming no periodic parameters for simplicity

  call nestRun( &
    nest_IS=.true., nest_mmodal=.false., nest_ceff=.false., &
    nest_nlive=Nlive, nest_tol=0.1d0, nest_ef=0.8d0, &
    nest_ndims=Ndims, nest_totPar=TotPar, nest_nCdims=Ndims, maxClst=10, &
    nest_updInt=100, nest_Ztol=-1.0d90, nest_root="my_chain_", &
    seed=-1, nest_pWrap=pWrap, nest_fb=.true., nest_resume=.false., &
    nest_outfile=.true., initMPI=.true., nest_logZero=-1.0d100, nest_maxIter=0, &
    loglike=my_loglikelihood, dumper=my_dumper, context=0 &
  )

contains

  subroutine my_loglikelihood(Cube, n_dim, nPar, lnew, context_pass)
    integer, intent(in) :: n_dim, nPar, context_pass
    double precision, intent(inout) :: Cube(nPar)
    double precision, intent(out) :: lnew
    ! ... calculate lnew based on Cube ...
    lnew = -0.5d0 * sum(Cube(1:n_dim)**2) ! Example: Gaussian
  end subroutine my_loglikelihood

  subroutine my_dumper(nSamples, nlive_d, nPar_d, physLive, posterior, paramConstr, &
                       maxLogLike, logZ, INSlogZ, logZerr, context_pass)
    integer, intent(in) :: nSamples, nlive_d, nPar_d, context_pass
    double precision, pointer, intent(in) :: physLive(:,:), posterior(:,:), paramConstr(:)
    double precision, intent(in) :: maxLogLike, logZ, INSlogZ, logZerr
    ! ... custom output ...
    print *, "Dumper called: logZ = ", logZ
  end subroutine my_dumper

end program my_analysis
```

## Dependencies and Interactions

*   **`utils1.f90`**: Uses various mathematical and matrix utilities (e.g., `Diagonalize`, `LogSumExp`, ellipsoid functions, random point generation within ellipsoids).
*   **`kmeans_clstr.f90`**: Called by `clusteredNest` (specifically `Dinosaur` routine which uses `Kmeans3`) for k-means clustering if multimodal sampling is active.
*   **`xmeans_clstr.f90`**: Called by `clusteredNest` (specifically `Dinosaur` routine) for X-means clustering (which can automatically determine K) if multimodal sampling is active.
*   **`posterior.F90`**: The `pos_samp` routine from this module is called at the end of `clusteredNest` or if `updInt` triggers it, to process the collected samples and generate posterior statistics.
*   **`priors.f90`**: While not directly calling functions from `priors.f90`, the `loglike` function provided by the user is responsible for transforming the unit hypercube samples (`Cube` array, whose elements are initially in [0,1]) into physical parameter values using appropriate prior transformations. The `getrandom` subroutine generates points in the unit hypercube.
*   **MPI Library**: Used for parallel execution if MultiNest is compiled with MPI support. `nestRun` handles MPI initialization and finalization if requested.
*   **C Wrapper (`cwrapper.f90`)**: If MultiNest is called from C, `cwrapper.f90` acts as an intermediary, translating C function calls and pointers to be compatible with `nestRun`.
```
