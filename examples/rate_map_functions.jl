#= Count number of trajectories in state i: =#

function Ni_trajectories(i_to_j_times,j_to_i_times,M,t)
    sum(j_to_i_times .≤ t) - sum(i_to_j_times .≤ t) + M
end



#= Construct the rate map: =#

function rate_map(b_to_u_forces,u_to_b_forces,dF,F_offset,M,N_bins)
    b_to_u_times = b_to_u_forces ./ dF
    u_to_b_times = (u_to_b_forces .+ F_offset) ./ dF
    histogram = fit(Histogram, b_to_u_forces, closed=:left, nbins=N_bins)
    F_0 = histogram.edges[1][1]
    ΔF = histogram.edges[1][2] - histogram.edges[1][1]
    C = histogram.weights
    N = length(C)
    b_to_u_rate_map = Array{Array{Float64,2},1}(undef, 0)
    for n = 1 : N
        t = histogram.edges[1][n+1]/dF
        N_b = Ni_trajectories(b_to_u_times,u_to_b_times,M,t)
        if C[n] > 0 && N_b > 0
            append!(b_to_u_rate_map,[[ (F_0 + (n-1/2)*ΔF) C[n]*dF/ΔF/N_b ]])
        end
    end
    histogram = fit(Histogram, u_to_b_forces, closed=:left, nbins=N_bins)
    F_0 = histogram.edges[1][1]
    ΔF = histogram.edges[1][2] - histogram.edges[1][1]
    C = histogram.weights
    N = length(C)
    u_to_b_rate_map = Array{Array{Float64,2},1}(undef, 0)
    for n = 1 : N
        t = (histogram.edges[1][n+1] + F_offset)/dF
        N_u = Ni_trajectories(u_to_b_times,b_to_u_times,0,t)
        if C[n] > 0 && N_u > 0
            append!(u_to_b_rate_map,[[ (F_0 + (n-1/2)*ΔF) C[n]*dF/ΔF/N_u ]])
        end
    end
    return vcat(b_to_u_rate_map...), vcat(u_to_b_rate_map...)
end