# RebindingKineticsMLE

This package provides a robust framework to analyze pulling traces from single-molecule force spectroscopy experiments involving multiple unbinding-rebinding events.  It makes use of numerically efficient maximum and asymptotically exact likelihood estimators for the parameters of the force-dependent unbinding and rebinding Bell-like rates.  

For more details on the theoretical framework, please refer to the associated preprint:
> J. T. Bullerjahn and G. Hummer, "Rebinding kinetics from single-molecule force spectroscopy experiments close to equilibrium", *arXiv:2205.05991* (2022). 
https://doi.org/10.48550/arXiv.2205.05991

Please cite the reference above if you use `RebindingKineticsMLE` to analyze your data.  



## Installation

The package is written in the open-source programming language [Julia](https://github.com/JuliaLang/julia), which can be downloaded from their [webpage](https://julialang.org/downloads/#download_julia).  

Currently, the package is not in a registry.  It must therefore be added by specifying a URL to the repository:
```julia
using Pkg; Pkg.add(url="https://github.com/bio-phys/RebindingKineticsMLE")
```
Users of older software versions may need to wrap the contents of the brackets with `PackageSpec()`.  



## Usage

### Importing data



### Parameter and error estimation

The Bell rate has two parameters, a rate and a length scale.  We therefore have four parameters, `Δx_off`, `k_off`, `Δx_on` and `k_on`, which can be estimated using the `MLE_estimator` function:
```julia
estimates, errors = MLE_estimator(data,Δt)
```
Here, `Δt` denotes the constant time step between two subsequent measurements.  
