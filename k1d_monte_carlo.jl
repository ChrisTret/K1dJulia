module K1dMonteCarlo

include("k1d_functionsMod.jl")

using .K1dFun, Distributions

export monte_carlo_sim, monte_carlo_test, extract_mcecb

"""
    monte_carlo_sim(data::Dict{String, Vector{Int64}}, key_pair::Tuple{String, String}, 
                    T::Vector{Int64}, m::Int64 = 1000) -> Matrix{Float64}

Performs a Monte Carlo simulation to generate random samples and compute univariate or bivariate statistics based on the provided data. The simulation can be used to assess the statistical properties of observed data by comparing it to randomly generated samples.

# Arguments
- `data::Dict{String, Vector{Int64}}`: A dictionary where each key is a string representing a dataset, and each value is a vector of integer data points.
- `key_pair::Tuple{String, String}`: A tuple containing two keys that specify which datasets in `data` to use for the simulation.
- `T::Vector{Int64}`: A vector of integer values representing the points at which the statistics will be computed.
- `m::Int64 = 1000`: The number of Monte Carlo iterations to perform (default is 1000).

# Returns
- `Matrix{Float64}`: A matrix where each column represents the results of one Monte Carlo simulation iteration. The number of rows corresponds to the length of `T`.
"""
function monte_carlo_sim(data::Dict{String, Vector{Int64}}, key_pair::Tuple{String, String}, T::Vector{Int64}, m::Int64 = 1000)
    key1 = key_pair[1]
    key2 = key_pair[2]
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
        area = maximum([maximum(X), maximum(Y)]) - minimum([minimum(X), minimum(Y)]) + 1    
        lambda_hat_x = Nx / area
        lambda_hat_y = Ny / area
        mc_results[:,1] = k1d_biv(X,Y,T)
        for i in 2:m
            mc_sample_x = round.(Int, cumsum(rand(Exponential(1/lambda_hat_x), Nx)))
            mc_sample_y = round.(Int, cumsum(rand(Exponential(1/lambda_hat_y), Ny)))
            mc_results[:,i] = k1d_biv(mc_sample_x, mc_sample_y, T)
        end
    end

    return mc_results
end


"""
    monte_carlo_sim(data::Dict{String, Dict{String, Dict{Int, Vector{Int64}}}}, 
                                             key_pair::Tuple{String, String}, T::Vector{Int64}, 
                                             m::Int64 = 1000) -> Dict{Tuple{String, String, String, Int}, Matrix{Float64}}

Performs a Monte Carlo simulation to generate random samples and compute univariate or bivariate statistics based on the provided genomic data separated by chromosome and region. The simulation assesses the statistical properties of observed data by comparing it to randomly generated samples.

# Arguments
- `data::Dict{String, Dict{String, Dict{Int, Vector{Int64}}}}`: A nested dictionary where each key is a group identifier, each sub-key is a chromosome identifier, and each sub-sub-key is a region index mapped to a vector of genomic positions.
- `key_pair::Tuple{String, String}`: A tuple containing two group names that specify which datasets in `data` to use for the simulation.
- `T::Vector{Int64}`: A vector of integer values representing the points at which the statistics will be computed.
- `m::Int64 = 1000`: The number of Monte Carlo iterations to perform (default is 1000).

# Returns
- `Dict{Tuple{String, String, String, Int}, Matrix{Float64}}`: A dictionary where each key is a tuple containing two group names, a chromosome identifier, and a region index. Each value is a matrix where each column represents the results of one Monte Carlo simulation iteration, and the number of rows corresponds to the length of `T`.
"""
function monte_carlo_sim(
    data::Dict{String, Dict{String, Dict{Int, Vector{Int64}}}}, 
    key_pair::Tuple{String, String}, 
    T::Vector{Int64}, 
    m::Int64 = 1000)::Dict{Tuple{String, String, String, Int}, Matrix{Float64}}
    key1 = key_pair[1]
    key2 = key_pair[2]
    mc_results_dict = Dict{Tuple{String, String, String, Int}, Matrix{Float64}}()

    # Loop through chromosomes
    for chrom in intersect(keys(data[key1]), keys(data[key2]))
        # Loop through regions
        for region in intersect(keys(data[key1][chrom]), keys(data[key2][chrom]))
            X = data[key1][chrom][region]
            Y = data[key2][chrom][region]
            mc_results = zeros(length(T), m)

            if key1 == key2
                N = length(X)
                lambda_hat = N / (X[N] - X[1] + 1)
                mc_results[:, 1] = k1d_univ(X, T)
                for i in 2:m
                    mc_sample = round.(Int, cumsum(rand(Exponential(lambda_hat), N)))
                    mc_results[:, i] = k1d_univ(mc_sample, T)
                end
            else
                Nx = length(X)
                Ny = length(Y)
                lambda_hat_x = Nx / (X[Nx] - X[1] + 1)
                lambda_hat_y = Ny / (Y[Ny] - Y[1] + 1)
                mc_results[:, 1] = k1d_biv(X, Y, T)
                for i in 2:m
                    mc_sample_x = round.(Int, cumsum(rand(Exponential(1 / lambda_hat_x), Nx)))
                    mc_sample_y = round.(Int, cumsum(rand(Exponential(1 / lambda_hat_y), Ny)))
                    mc_results[:, i] = k1d_biv(mc_sample_x, mc_sample_y, T)
                end
            end

            mc_results_dict[(key1, key2, chrom, region)] = mc_results
        end
    end

    return mc_results_dict
end



"""
    extract_mcecb(sim_vals::Matrix{Float64}, alpha::Float64 = 0.05) -> Matrix{Float64}

Calculates the Monte Carlo error confidence bounds (MCECB) for each row of the simulation values matrix.

# Arguments
- `sim_vals::Matrix{Float64}`: A matrix of simulation values where each row represents a different simulation scenario and each column (except the first) represents a different simulation iteration.
- `alpha::Float64 = 0.05`: The significance level used to calculate the lower and upper percentiles for the confidence bounds (default is 0.05).

# Returns
- `Matrix{Float64}`: A matrix with two columns where the first column contains the lower percentiles and the second column contains the upper percentiles for each row of the input matrix.
"""
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


"""
    monte_carlo_test(sim_vals::Matrix{Float64}, alpha::Float64 = 0.05) -> Vector{Int}

Performs a Monte Carlo test to identify rows in the simulation values matrix where the first column value falls outside the specified confidence bounds.

# Arguments
- `sim_vals::Matrix{Float64}`: A matrix of simulation values where each row represents a different simulation scenario. The first column contains observed values, and the remaining columns contain simulation results.
- `alpha::Float64 = 0.05`: The significance level used to calculate the lower and upper percentiles for determining the confidence bounds (default is 0.05).

# Returns
- `Vector{Int}`: A vector of row indices where the observed value (first column) is outside the calculated confidence bounds based on the simulation results.
"""
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


# # Create the dictionary
# bias_dict = Dict{String, Vector{Int64}}()

# # Populate the dictionary with the specified keys and values
# bias_dict["big_lambda"] = vcat(1:999, 50_000)
# bias_dict["medium_lambda"] = vcat(1:999, 500_000)
# bias_dict["small_lambda"] = vcat(1:999, 5_000_000)

end # end of module