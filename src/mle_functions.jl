#= Sort the data into the possible four i –> j transitions: =#

function sort_data(data#=::AbstractArray{Array{Float64,2},1}=#)
    M = size(data,1) # Number of time series
    N = size.(data,1) # Array of lengths of time series
    b_to_b = Array{Float64,1}(undef,0)
    u_to_b = Array{Float64,1}(undef,0)
    b_to_u = Array{Float64,1}(undef,0)
    u_to_u = Array{Float64,1}(undef,0)
    for m = 1 : M
        tmp_data = data[m]
        for n = 2 : N[m]
            if tmp_data[n,2] == 0.0 # bound state
                if tmp_data[n-1,2] == 0.0 # bound state
                    push!(b_to_b,tmp_data[n-1,1])
                else
                    push!(u_to_b,tmp_data[n-1,1])
                end
            else
                if tmp_data[n-1,2] == 0.0 # bound state
                    push!(b_to_u,tmp_data[n-1,1])
                else
                    push!(u_to_u,tmp_data[n-1,1])
                end
            end
        end
    end
    return b_to_b, u_to_b, b_to_u, u_to_u
end



#= Maximum likelihood estimator for the instantaneous rate k_off^0: =#

function k_bu0(β,Δt,Δxb,N_bu,N_bb,b_to_b)
#    N_bb = length(b_to_b)
    tmp = 0.0
    for n = 1 : N_bb
        tmp += exp(β*Δxb*b_to_b[n])
    end
    k_bu0 = N_bu / (Δt*tmp)
    return k_bu0
end

#= Reduced negative-log-likelihood for Δx_off: =#

function Δxb_likelihood(β,Δxb,N_bu,sum_b_to_u,N_bb,b_to_b)
#    N_bb = length(b_to_b)
    tmp = 0.0
    for n = 1 : N_bb
        tmp += exp(β*Δxb*b_to_b[n])
    end
    L = N_bu*log(tmp) - β*Δxb*sum_b_to_u
    return L / N_bb
end

#= Maximum likelihood estimator for the instantaneous rate k_on^0: =#

function k_ub0(β,Δt,Δxu,N_ub,N_uu,u_to_u)
#    N_uu = length(u_to_u)
    tmp = 0.0
    for n = 1 : N_uu
        tmp += exp(-β*Δxu*u_to_u[n])
    end
    k_ub0 = N_ub / (Δt*tmp)
    return k_ub0
end

#= Reduced negative-log-likelihood for Δx_on: =#

function Δxu_likelihood(β,Δxu,N_ub,sum_u_to_b,N_uu,u_to_u)
#    N_uu = length(u_to_u)
    tmp = 0.0
    for n = 1 : N_uu
        tmp += exp(-β*Δxu*u_to_u[n])
    end
    L = N_ub*log(tmp) + β*Δxu*sum_u_to_b
    return L / N_uu
end



#= Evaluate the uncertainty of the estimates: =#

function bu_errors(β,Δt,Δxb,k_bu0,N_bu,N_bb,b_to_b,δF)
    I_bu_11 = 0.0
    I_bu_12 = 0.0
    for n = 1 : N_bb
        tmp = exp(β*Δxb*b_to_b[n])
        I_bu_11 += tmp*(δF^2 + (b_to_b[n] + β*Δxb*δF^2)^2)
        I_bu_12 += tmp*(b_to_b[n] + β*Δxb*δF^2)
    end
    I_bu_11 = β^2*k_bu0*Δt*exp(0.5*(β*Δxb*δF)^2)*I_bu_11
    I_bu_12 = β*Δt*exp(0.5*(β*Δxb*δF)^2)*I_bu_12
    I_bu_22 = N_bu/k_bu0^2
    det_bu = I_bu_11*I_bu_22 - I_bu_12^2
    δΔxb = sqrt(I_bu_22/det_bu)
    δk_bu0 = sqrt(I_bu_11/det_bu)
    return δΔxb, δk_bu0
end

function ub_errors(β,Δt,Δxu,k_ub0,N_ub,N_uu,u_to_u,δF)
    I_ub_11 = 0.0
    I_ub_12 = 0.0
    for n = 1 : N_uu
        tmp = exp(-β*Δxu*u_to_u[n])
        I_ub_11 += tmp*(δF^2 + (u_to_u[n] - β*Δxu*δF^2)^2)
        I_ub_12 += tmp*(u_to_u[n] - β*Δxu*δF^2)
    end
    I_ub_11 = β^2*k_ub0*Δt*exp(0.5*(β*Δxu*δF)^2)*I_ub_11
    I_ub_12 = -β*Δt*exp(0.5*(β*Δxu*δF)^2)*I_ub_12
    I_ub_22 = N_ub/k_ub0^2
    det_ub = I_ub_11*I_ub_22 - I_ub_12^2
    δΔxu = sqrt(I_ub_22/det_ub)
    δk_ub0 = sqrt(I_ub_11/det_ub)
    return δΔxu, δk_ub0
end



#= Compute MLE estimates for the parameters Δx_off, k_off^0, Δx_on, and k_on^0: =#

function MLE_estimator(data,Δt,β=1/4,δF=0.0,interval=[0.0,10.0])
    b_to_b, u_to_b, b_to_u, u_to_u = sort_data(data)
    N_bb = length(b_to_b)
    N_ub = length(u_to_b)
    N_bu = length(b_to_u)
    N_uu = length(u_to_u)
    sum_b_to_u = sum(b_to_u)
    sum_u_to_b = sum(u_to_b)
    # In this part, we estimate the parameters:
    res_1 = optimize(Δxb->Δxb_likelihood(β,Δxb,N_bu,sum_b_to_u,N_bb,b_to_b), interval[1], interval[2], Brent())
    Δxb_mle = res_1.minimizer
    k_bu0_mle = k_bu0(β,Δt,Δxb_mle,N_bu,N_bb,b_to_b)
    res_2 = optimize(Δxu->Δxu_likelihood(β,Δxu,N_ub,sum_u_to_b,N_uu,u_to_u), interval[1], interval[2], Brent())
    Δxu_mle = res_2.minimizer
    k_ub0_mle = k_ub0(β,Δt,Δxu_mle,N_ub,N_uu,u_to_u)
    # In this part, we calculate the associated error:
    δΔxb_mle, δk_bu0_mle = bu_errors(β,Δt,Δxb_mle,k_bu0_mle,N_bu,N_bb,b_to_b,δF)
    δΔxu_mle, δk_ub0_mle = ub_errors(β,Δt,Δxu_mle,k_ub0_mle,N_ub,N_uu,u_to_u,δF)
    return [Δxb_mle, Δxu_mle, k_bu0_mle, k_ub0_mle], [ δΔxb_mle, δΔxu_mle, δk_bu0_mle, δk_ub0_mle ]
end
