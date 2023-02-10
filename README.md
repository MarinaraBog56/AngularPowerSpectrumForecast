# AngularPowerSpectrumForecast
A collection of codes used to generate Fisher forecasts for 3x2pt surveys.

## Requirements
The jupyter notebook codes require CAMB 1.4.0 to run, and a valid FILEPATH stated at the top of each notebook.
The C++ code requries GSL-2.7 to compile and run, and a valid FILEPATH stated in the stdafx.h file.
Computing the 4MOST galaxy selection functions requires cleaned data from the 4MOST Facility Simulations.

## Usage
Codes should be ran in order:
-InputGeneration
-AngularSpectraGeneration
-FisherMatrixGeneration
-GraphGeneration

InputGeneration contains a jupyter notebook for generating the input files required to compute angular power spectra.
Additionally, includes:
- Functions for checking galaxy selection functions as well as their corresponding kernel integrands.
- Functions for computing transfer functions.

AngularSpectraGeneration contains C++ code for calculating all combinations of angular power spectra required for the Fisher forecast. Code allows for multi-threading, making it ideal for usage in a HPC.
- Changing the number of threads is a command line argument.
- It is recommended to split this code up by parameters, and run as multiple jobs so as to further quicken the computation.
- Computing 10 parameter deviations + mean exactly from 0<ell<50, and Limber approximation for ell>50 with 8 threads takes ~5 hours.

FisherMatrixGeneration contains a jupyter notebook for computing Fisher matrices given instrument models and angular power spectra.
- Also includes functions which return signal-to-noise plots and correlation matrix plots.
- See the code for examples of this.

GraphGeneration contains a jupyter notebook for graphing Fisher matrices into corner contour plots.
- Additionally includes functions for retrieving and graphing the Figure of Merit.

## Other notes
C++ code solves radial integrals through ODE solvers, and then splining the final integrand. This can produce some noise at high angular multipole moments, unless you use fine sampling of k-space. Possible workarounds for this include implementing FFTLog.
The parameter file is saved as a CSV. This occasionally introduces a bug, requiring you to open the CSV file and save it to resolve.

## Future planned updates
- Remove the CMB Lensing kernels from computation.
- Fold in galaxyLens and galaxyClustering kernels into a single computation.
- Long term: implement FFTLog, wrap C code with wrapper, implement MCMC routine.

## License
[MIT](https://choosealicense.com/licenses/mit/)
