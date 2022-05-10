# Code for estimating the parameters of Bell-like unbinding and rebinding rates
# version 1.0 (10/05/2022)
# Jakob Tómas Bullerjahn (jabuller@biophys.mpg.de)
# Gerhard Hummer

# Please read and cite the associated publication: 
# J. T. Bullerjahn and G. Hummer, "Rebinding kinetics from single-molecule force spectroscopy experiments close to equilibrium", ...



module RebindingKineticsMLE

export
        # MLE:
        MLE_estimator,
        sort_data

using Optim

include("mle_functions.jl")

end
