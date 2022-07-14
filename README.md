# RebindingKineticsMLE

This package provides a robust framework to analyze pulling traces from single-molecule force spectroscopy experiments involving multiple unbinding-rebinding events.  It makes use of numerically efficient and asymptotically exact maximum likelihood estimators for the parameters of the force-dependent unbinding and rebinding Bell-like rates.  

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
Users of older versions of Julia may need to wrap the contents of the brackets with `PackageSpec()`.  



## Usage

A detailed example on how to apply the code and compare its output to rate maps can be found in [the examples directory](examples).  



### Importing data

Every pulling trace (`Array{Float64,2}`) is defined by two columns.  The first column lists the applied force (in *pN*) at each time step, whereas the second column defines the associated state the system occupies (`0.0` for bound state, `1.0` for unbound state).  A collection of pulling traces should therefore be of the type `Array{Array{Float64,2},1}`.  

Assuming that each pulling trace is saved in a separate file, stored in the directory `pulling_traces`, we can use something like the following code snippet to import the data:
```julia
using DelimitedFiles

function read_data(dir)
    files = readdir(dir)
    N = length(files)
    data = Array{Array{Float64,2},1}(undef, N)
    for n = 1 : N
        data[n] = readdlm(string(dir,files[n]))
    end
    return data
end

data = read_data("./pulling_traces/")
```
Note that the pulling traces do not have to be of equal length or have been generated using the same force protocol.  The only requirement is that the time step `Δt` between two subsequent force measurements has to be constant and the same for all pulling traces.  



### Parameter and error estimation

The Bell rate has two parameters, a rate and a length scale.  We therefore have four parameters, `Δx_off`, `k_off`, `Δx_on` and `k_on` (lengths in *nm*, rates in *1/s*), which can be estimated using the `MLE_estimator` function:
```julia
estimates, errors = MLE_estimator(data,Δt)
```
for a constant time step `Δt` between two subsequent measurements.  `MLE_estimator` has a few optional arguments: 
```julia
MLE_estimator(data,Δt,β=1/4,δF=0.0,interval=[0.0,10.0])
```
where `β=1/(k_B*T)` [in *1/(pN nM)*] denotes the inverse thermal energy scale with the Boltzmann constant `k_B` and temperature `T`, `δF` characterizes the experimentally determined mean-squared fluctuations of the applied force, which we assume to be constant, and `interval` sets the search range used by the [optimizer](https://github.com/JuliaNLSolvers/Optim.jl) to determine `Δx_off` and `Δx_on` (both in `nm`).  
