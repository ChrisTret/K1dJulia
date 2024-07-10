module K1dMonteCarlo

include("k1d_functionsMod.jl")

using .K1dFun, Distributions

export monte_carlo_sim

function monte_carlo_sim(data::Dict{String,Vector{Int64}}, keys::Tuple{String, String}, T::Vector{Int64}, m::Int64 = 100)
    key1 = keys[1]
    key2 = keys[2]
    X = data[key1]
    Y = data[key2]
    mc_results = zeros(length(T),m)

    if key1 == key2
        N = length(X)
        lambda_hat = N / (X[N] - X[1] + 1)
        mc_results[:,1] = k1d_univ(X,T)
        for i in 2:m
            mc_sample = round.(Int, cumsum(rand(Exponential(1/lambda_hat), N)))
            mc_results[:,i] = k1d_univ(mc_sample, T)
        end
    else
        Nx = length(X)
        Ny = length(Y)
        lambda_hat_x = Nx / (X[Nx] - X[1] + 1)
        lambda_hat_y = Ny / (Y[Ny] - Y[1] + 1)
        mc_results[:,1] = k1d_biv(X,Y,T)
        for i in 2:3
            mc_sample_x = round.(Int, cumsum(rand(Exponential(1/lambda_hat_x), Nx)))
            mc_sample_y = round.(Int, cumsum(rand(Exponential(1/lambda_hat_y), Ny)))
            mc_results[:,i] = k1d_biv(mc_sample_x, mc_sample_y, T)
        end
    end

    return mc_results
end


function monte_carlo_test(sim_vals, alpha::Float64 = 0.05)


end

end # end of module