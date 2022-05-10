const β = 1/4 # (pN*nm)^(-1)

#= Rate expressions and related functions: =#

k_on(Δx_on,k_on0,F,F_0) = k_on0 * exp(-β*Δx_on*(F-F_0))
k_off(Δx_off,k_off0,F) = k_off0 * exp(β*Δx_off*F)



#= Gillespie algorithm: =#

function t_on(Δx_on,k_on0,dF,F_0,t_prev,R)
    tmp = k_on(Δx_on,k_on0,dF*t_prev,F_0) + β*Δx_on*dF*log(R)
    if tmp > 0.0
        return log(k_on0/tmp)/(β*Δx_on*dF) + F_0/dF
    else
        return Inf
    end
end

function t_off(Δx_off,k_off0,dF,t_prev,R)
    tmp = k_off(Δx_off,k_off0,dF*t_prev) - β*Δx_off*dF*log(R)
    return log(tmp / k_off0) / (β*Δx_off*dF)
end

function generate_transition_times(dF,F_0,parameters)
    Δx_off = parameters[1]
    Δx_on = parameters[2]
    k_off0 = parameters[3]
    k_on0 = parameters[4]
    τ_off = Array{Float64,1}(undef, 0)
    τ_on = Array{Float64,1}(undef, 0)
    s = 0.0
    t_prev = 0.0
    while t_prev < Inf
        if s == 0.0
            R = rand()
            t_prev = t_off(Δx_off,k_off0,dF,t_prev,R)
            append!(τ_off,t_prev)
            s = 1.0
        else
            R = rand()
            t_prev = t_on(Δx_on,k_on0,dF,F_0,t_prev,R)
            if t_prev < Inf
                append!(τ_on,t_prev)
                s = 0.0
            end
        end
    end
    return τ_off, τ_on
end

function make_τ(M,dF,F_0,parameters)
    τ_off = Array{Array{Float64,1},1}(undef, M)
    τ_on = Array{Array{Float64,1},1}(undef, M)
    for m = 1 : M
        τ_off[m], τ_on[m] = generate_transition_times(dF,F_0,parameters)
    end
    return τ_off, τ_on
end

function N_correction(M_off, M_on) # This does not account for multiple datasets M
    M = size(M_off,1)
    N = length.(M_on)
    N_off = [ Array{Int64,1}(undef,0) for m = 1 : M ]
    N_on = [ Array{Int64,1}(undef,0) for m = 1 : M ]
    for m = 1 : M
        append!(N_off[m],M_off[m][1])
        for i = 1 : N[m]
            if M_off[m][i+1] > M_on[m][i]
                append!(N_off[m],M_off[m][i+1])
                append!(N_on[m],M_on[m][i])
            end
        end
    end
    return N_off, N_on
end

function make_data(M,dF,F_0,Δt,parameters)
    data = Array{Array{Float64,2},1}(undef, M)
    τ_off, τ_on = make_τ(M,dF,F_0,parameters)
    N_off_tmp = [ round.(Int64,τ_off[m]/Δt) .+ 1 for m = 1 : M ]
    N_on_tmp = [ round.(Int64,τ_on[m]/Δt) .+ 1 for m = 1 : M ]
    N_off, N_on = N_correction(N_off_tmp, N_on_tmp)
    for m = 1 : M
        N = length(N_on[m])
        tmp_data = zeros(5*N_off[m][1] + N_off[m][end],2)
        for i = 1 : N_off[m][1]
            tmp_data[i,1] = dF*(i-1)*Δt
            tmp_data[i,2] = 0.0
        end
        for n = 1 : N
            for i = N_off[m][n]+1 : N_on[m][n]
                tmp_data[i,1] = dF*(i-1)*Δt - F_0
                tmp_data[i,2] = 1.0
            end
            for i = N_on[m][n]+1 : N_off[m][n+1]
                tmp_data[i,1] = dF*(i-1)*Δt
                tmp_data[i,2] = 0.0
            end
        end
        for i = N_off[m][end]+1 : 5*N_off[m][1] + N_off[m][end]
            tmp_data[i,1] = dF*(i-1)*Δt - F_0
            tmp_data[i,2] = 1.0
        end
        data[m] = tmp_data
    end
    return data
end
