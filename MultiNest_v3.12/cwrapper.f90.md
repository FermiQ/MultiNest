# Documentation for `cwrapper.f90`

## Overview

This file provides a C language wrapper for the core Fortran routines in MultiNest. It allows users to call MultiNest's nested sampling algorithm (`nestRun`) from C or C++ applications. The interfacing is achieved using the `iso_c_binding` standard Fortran module.

## Key Components

### `loglike_f`
*   **Type**: Subroutine
*   **Description**: This is a Fortran wrapper for the C log-likelihood function provided by the user. It translates the C function pointer and arguments into a Fortran-callable interface.
*   **Parameters**:
    *   `Cube`: (Input) Array of parameter values.
    *   `n_dim`: (Input) Number of dimensions (parameters).
    *   `nPar`: (Input) Total number of parameters.
    *   `lnew`: (Output) The calculated log-likelihood value.
    *   `context_pass`: (Input) A context pointer (integer) passed from C.

### `dumper_f`
*   **Type**: Subroutine
*   **Description**: This is a Fortran wrapper for the C dumper function provided by the user. The dumper function is called by MultiNest at regular intervals to output live points or other information.
*   **Parameters**:
    *   `nSamples`: (Input) Number of samples in `posterior`.
    *   `nlive`: (Input) Number of live points.
    *   `nPar`: (Input) Total number of parameters.
    *   `physLive`: (Input, Pointer) Array of physical live points.
    *   `posterior`: (Input, Pointer) Array of posterior samples.
    *   `paramConstr`: (Input, Pointer) Array of parameter constraints.
    *   `maxLogLike`: (Input) Maximum log-likelihood value.
    *   `logZ`: (Input) Current estimate of the log-evidence.
    *   `INSlogZ`: (Input) Current estimate of the importance nested sampling log-evidence.
    *   `logZerr`: (Input) Current estimate of the error in log-evidence.
    *   `context_pass`: (Input) A context pointer (integer) passed from C.

### `run`
*   **Type**: Subroutine (`bind(c)`)
*   **Description**: This is the main C-callable subroutine that initiates the MultiNest nested sampling algorithm. It takes various control parameters, the C function pointers for log-likelihood and dumper, and a context pointer.
*   **Parameters (主なもの)**:
    *   `nest_IS`: (Input, `logical(c_bool)`) Importance Sampling on/off.
    *   `nest_mmodal`: (Input, `logical(c_bool)`) Multimodal sampling on/off.
    *   `nest_ceff`: (Input, `logical(c_bool)`) Constant efficiency mode on/off.
    *   `nest_nlive`: (Input, `integer(c_int)`) Number of live points.
    *   `nest_tol`: (Input, `real(c_double)`) Tolerance for evidence calculation.
    *   `nest_ef`: (Input, `real(c_double)`) Efficiency factor.
    *   `nest_ndims`: (Input, `integer(c_int)`) Number of dimensions (parameters being varied).
    *   `nest_totPar`: (Input, `integer(c_int)`) Total number of parameters (including derived).
    *   `nest_nCdims`: (Input, `integer(c_int)`) Number of parameters to be clustered.
    *   `nest_updInt`: (Input, `integer(c_int)`) Update interval for output files.
    *   `nest_root`: (Input, `character(kind=c_char,len=1) dimension(1)`) Root name for output files.
    *   `seed`: (Input, `integer(c_int)`) Random seed.
    *   `loglike`: (Input, `type(c_funptr)`) C function pointer to the log-likelihood function.
    *   `dumper`: (Input, `type(c_funptr)`) C function pointer to the dumper function.
    *   `context`: (Input, `type(c_ptr)`) Context pointer passed to `loglike` and `dumper`.

## Important Variables/Constants

*   **`theloglike`**: Module variable of `type(c_funptr)` storing the C log-likelihood function pointer.
*   **`thedumper`**: Module variable of `type(c_funptr)` storing the C dumper function pointer.

## Usage Examples

To use MultiNest from C, a user would typically:
1.  Define a log-likelihood function in C with the prototype:
    `double loglike(double *Cube, int n_dim, int nPar, void *context);`
2.  Optionally, define a dumper function in C with the prototype:
    `void dumper(int nSamples, int nlive, int nPar, double *physLive, double *posterior, double *paramConstr, double maxLogLike, double logZ, double INSlogZ, double logZerr, void *context);`
3.  Call the `run` subroutine from C, passing the function pointers and other necessary parameters.

```c
// Example C-style call (conceptual)
#include "multinest.h" // Assuming multinest.h would declare the 'run' function

// ... (define loglikelihood_c and dumper_c) ...

void *my_context = NULL; // Or some custom context data

run(true, true, false, 1000, 0.5, 0.8, Ndim, Npar, Npar_cluster, 100,
    1000, -1e90, "chains/test-", -1, pWrap_c,
    true, false, true, true, -1e100, 0,
    loglikelihood_c, dumper_c, my_context);
```

## Dependencies and Interactions

*   **`iso_c_binding`**: Standard Fortran module used for C interoperability.
*   **`Nested` module**: Specifically, the `nestRun` subroutine from the `Nested` module (likely in `nested.F90`) is called by the `run` subroutine in this wrapper.
```
