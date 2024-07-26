module K1dMonteCarlo

include("k1d_functionsMod.jl")

using .K1dFun, Distributions

export monte_carlo_sim, monte_carlo_test, extract_mcecb

function monte_carlo_sim(data::Dict{String,Vector{Int64}}, keys::Tuple{String, String}, T::Vector{Int64}, m::Int64 = 1000)
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
        for i in 2:m
            mc_sample_x = round.(Int, cumsum(rand(Exponential(1/lambda_hat_x), Nx)))
            mc_sample_y = round.(Int, cumsum(rand(Exponential(1/lambda_hat_y), Ny)))
            mc_results[:,i] = k1d_biv(mc_sample_x, mc_sample_y, T)
        end
    end

    return mc_results
end

function extract_mcecb(sim_vals::Matrix{Float64}, alpha::Float64 = 0.05)
    # Get the number of rows and columns
    num_rows = size(sim_vals, 1)

    # Initialize matrices to store the lower and upper percentiles
    lower_percentiles = zeros(num_rows)
    upper_percentiles = zeros(num_rows)

    # Iterate over each row
    for i in 1:num_rows
        # Extract the i-th row and sort it (excluding the first column)
        sorted_row = sort(sim_vals[i, 2:end])

        # Calculate the percentiles
        lower_percentiles[i] = quantile(sorted_row, alpha/2)
        upper_percentiles[i] = quantile(sorted_row, 1-alpha/2)
    end

    # Combine the lower and upper percentiles into a matrix with two columns
    mcecb = hcat(lower_percentiles, upper_percentiles)

    return mcecb
end

function monte_carlo_test(sim_vals::Matrix{Float64}, alpha::Float64 = 0.05)
    # Get the number of rows
    num_rows = size(sim_vals, 1)

    # Initialize an empty array to store rows that meet the condition
    significant_rows = []

    # Iterate over each row
    for i in 1:num_rows
        # Extract the i-th row
        row = sim_vals[i, :]

        # Sort the row to find the percentiles
        sorted_row = sort(row[2:end])

        # Calculate the percentiles
        lower_percentile = quantile(sorted_row, alpha/2)
        upper_percentile = quantile(sorted_row, 1-alpha/2)

        # Check if the first column value is outside the percentile range
        if row[1] < lower_percentile || row[1] > upper_percentile
            push!(significant_rows, i)
        end
    end

    return significant_rows
end

end # end of module