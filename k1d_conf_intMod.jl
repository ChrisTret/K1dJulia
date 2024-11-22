module K1dConfInt

include("k1d_functionsMod.jl")

using .K1dFun, Statistics, Random, Base.Threads

export bootstrap_variance_k1d, bootstrap_variance_k1d_all_comparisons, approx_k1d_mean_variance_chromosome, approx_k1d_mean_variance_all_comparisons

"""
    bootstrap_variance_k1d(Data::Dict{String, Vector{Int64}}, keys::Tuple{String,String}, T::Vector{Int64}, n::Int, B::Int) -> Dict{Tuple{String, String}, Vector{Float64}}

Calculate the bootstrap variance of the estimated K function.

# Arguments
- `Data:: Dict{String, Vector{Float64}}`: A dictionary where each key is a name and the value is a vector of `genoCenter` values.
- `keys::Tuple{String, String}`: Tuple of keys to access the data
- `T::Vector{Int64}`: A vector of values for the diameters.
- `B::Int`: The number of bootstrap samples.

# Returns
- `Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary where each key is a tuple of names from `Data` and the value is the bootstrap variance vector.
"""
function bootstrap_variance_k1d(Data::Dict{String, Vector{Int64}}, keys::Tuple{String,String}, T::Vector{Int64}, B::Int)
    bootstrap_results = Dict{Tuple{String, String}, Vector{Float64}}()
    K_bootstrap = zeros(Float64, length(T), B)
    X = Data[keys[1]]
    Y = Data[keys[2]]

    if keys[1] == keys[2]
        X_num = length(X)
        @threads for b in 1:B
            sample_X = rand(X,X_num)
            K_bootstrap[:, b] = k1d_univ(sample_X, T)
        end
    else
        X_num = length(Data[keys[1]])
        Y_num = length(Data[keys[2]])
        @threads for b in 1:B
            sample_X = rand(X,X_num)
            sample_Y = rand(Y, Y_num)
            K_bootstrap[:, b] = k1d_biv(sample_X, sample_Y, T)
        end
    end

    variance_K = var(K_bootstrap, dims=2)
    bootstrap_results[(keys[1], keys[2])] = vec(variance_K)
    
    return bootstrap_results
end


"""
    bootstrap_variance_k1d_all_comparisons(Data::Dict{String, Vector{Int64}}, T::Vector{Int64}, n::Int, B::Int) -> Dict{Tuple{String, String}, Vector{Float64}}

Calculate the bootstrap variance of the estimated K function.

# Arguments
- `Data:: Dict{String, Vector{Float64}}`: A dictionary where each key is a name and the value is a vector of `genoCenter` values.
- `T::Vector{Int64}`: A vector of values for the diameters.
- `B::Int`: The number of bootstrap samples.

# Returns
- `Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary where each key is a tuple of names from `Data` and the value is the bootstrap variance vector.
"""
function bootstrap_variance_k1d_all_comparisons(Data::Dict{String, Vector{Int64}}, T::Vector{Int64}, B::Int)::Dict{Tuple{String, String}, Vector{Float64}}
    bootstrap_results = Dict{Tuple{String, String}, Vector{Float64}}()
    dict_keys = collect(keys(Data))
    
    # Initialize the parallel processing
    for key1 in dict_keys
        for key2 in dict_keys
            bootstrap_results[key1,key2] = bootstrap_variance_k1d(Data,(key1,key2), T, B)
        end
    end
    
    return bootstrap_results
end


"""
    approx_k1d_mean_variance_chromosome(results::Dict{Tuple{String, String, String}, Vector{Float64}}, 
                                        k_means::Dict{Tuple{String, String}, Vector{Float64}}, 
                                        key_pair::Tuple{String, String}) -> Dict{Tuple{String, String}, Vector{Float64}}

Calculates the approximate variance of the k1d mean for a specific key pair across different chromosomes based on provided results.

# Arguments
- `results::Dict{Tuple{String, String, String}, Vector{Float64}}`: A dictionary where each key is a tuple containing three strings (e.g., representing group, chromosome, and some identifier), and each value is a vector of float results.
- `k_means::Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary where each key is a tuple of two strings (e.g., representing group and another identifier), and each value is a vector representing the mean values of k1d.
- `key_pair::Tuple{String, String}`: A tuple containing two strings to specify which group and identifier in `k_means` to use for calculating variance.

# Returns
- `Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary where the key is the provided `key_pair` and the value is a vector of variances for the k1d mean.
"""
function approx_k1d_mean_variance_chromosome(results::Dict{Tuple{String, String, String}, Vector{Float64}}, k_means::Dict{Tuple{String, String}, Vector{Float64}}, key_pair::Tuple{String, String})
    len_T = length(k_means[key_pair])
    k_mean = k_means[key_pair]
    k_mean_var = Dict{Tuple{String, String}, Vector{Float64}}(key_pair => zeros(len_T))

    filtered_keys = filter(key -> key[1] == key_pair[1] && key[2] == key_pair[2], keys(results))
    R = length(filtered_keys)

    if R > 1
        for key in filtered_keys
            diff = results[key] .- k_mean
            k_mean_var[key_pair] .+= diff .* diff
        end

        scale_factor = 1 / (R * (R - 1))
        k_mean_var[key_pair] .*= scale_factor
    else
        k_mean_var[key_pair] .= NaN
        println("Not enough data points to calculate variance for key pair $(key_pair).")
    end

    return k_mean_var
end

"""
    approx_k1d_mean_variance_region(results::Dict{Tuple{String, String, String, Int}, Vector{Float64}}, 
                                    k_means::Dict{Tuple{String, String}, Vector{Float64}}, 
                                    key_pair::Tuple{String, String}) -> Dict{Tuple{String, String}, Vector{Float64}}

Calculates the approximate variance of the Kmean for a specific key pair across different regions based on the provided results.

# Arguments
- `results::Dict{Tuple{String, String, String, Int}, Vector{Float64}}`: 
  A dictionary where each key is a tuple containing two element names, a chromosome identifier, and a region identifier, and each value is a vector of K statistics for that region.
- `k_means::Dict{Tuple{String, String}, Vector{Float64}}`: 
  A dictionary where each key is a tuple of two element names, and each value is a vector representing the mean values of K across regions.
- `key_pair::Tuple{String, String}`: 
  A tuple containing two element names to specify which pair's variance to compute.

# Returns
- `Dict{Tuple{String, String}, Vector{Float64}}`: 
  A dictionary where the key is the provided `key_pair`, and the value is a vector of variances for the K mean across regions.
"""
function approx_k1d_mean_variance_region(results::Dict{Tuple{String, String, String, Int}, Vector{Float64}}, k_means::Dict{Tuple{String, String}, Vector{Float64}}, key_pair::Tuple{String, String})
    len_T = length(k_means[key_pair])
    k_mean = k_means[key_pair]
    k_mean_var = Dict{Tuple{String, String}, Vector{Float64}}(key_pair => zeros(len_T))

    # Filter keys matching the key_pair across all regions
    filtered_keys = filter(key -> key[1] == key_pair[1] && key[2] == key_pair[2], keys(results))
    R = length(filtered_keys)

    if R > 1
        for key in filtered_keys
            diff = results[key] .- k_mean
            k_mean_var[key_pair] .+= diff .* diff
        end

        scale_factor = 1 / (R * (R - 1))
        k_mean_var[key_pair] .*= scale_factor
    else
        k_mean_var[key_pair] .= NaN
        println("Not enough data points to calculate variance for key pair $(key_pair).")
    end

    return k_mean_var
end



"""
    approx_k1d_mean_variance_all_comparisons(results::Dict{Tuple{String, String, String}, Vector{Float64}}, 
                                             k_means::Dict{Tuple{String, String}, Vector{Float64}}) 
                                             -> Dict{Tuple{String, String}, Vector{Float64}}

Calculates the approximate variance of the k1d mean for all unique key pairs across different chromosomes based on the provided results.

# Arguments
- `results::Dict{Tuple{String, String, String}, Vector{Float64}}`: A dictionary where each key is a tuple containing three strings (e.g., representing group, chromosome, and some identifier), and each value is a vector of float results.
- `k_means::Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary where each key is a tuple of two strings (e.g., representing group and another identifier), and each value is a vector representing the mean values of k1d.

# Returns
- `Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary where each key is a tuple of two strings (representing a group and identifier) and each value is a vector of variances for the k1d mean calculated for that key pair.
"""
function approx_k1d_mean_variance_all_comparisons(results::Dict{Tuple{String, String, String}, Vector{Float64}}, k_means::Dict{Tuple{String, String}, Vector{Float64}})
    # Initialize the dictionary to store variances
    all_variances = Dict{Tuple{String, String}, Vector{Float64}}()

    # Iterate over all unique key pairs from the results dictionary
    for key in keys(results)
        key_pair = (key[1], key[2])
        
        # Check if the variance for this key pair has already been calculated
        if !haskey(all_variances, key_pair)
            variance = approx_k1d_mean_variance_chromosome(results, k_means, key_pair)
            all_variances[key_pair] = variance[key_pair]
        end
    end

    return all_variances
end


function approx_k1d_mean_variance_all_comparisons(
    results::Dict{Tuple{String, String, String, Int}, Vector{Float64}},
    k_means::Dict{Tuple{String, String}, Vector{Float64}}
)::Dict{Tuple{String, String}, Vector{Float64}}
    # Initialize the dictionary to store variances
    all_variances = Dict{Tuple{String, String}, Vector{Float64}}()

    # Iterate over all unique key pairs from the results dictionary
    for key in keys(results)
        key_pair = (key[1], key[2])
        
        # Check if the variance for this key pair has already been calculated
        if !haskey(all_variances, key_pair)
            variance = approx_k1d_mean_variance_region(results, k_means, key_pair)
            all_variances[key_pair] = variance[key_pair]
        end
    end

    return all_variances
end


end # end of module