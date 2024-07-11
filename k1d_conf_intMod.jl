module K1dConfInt

include("k1d_functionsMod.jl")

using .K1dFun, Statistics, Random, Base.Threads

export bootstrap_variance_k1d, bootstrap_variance_k1d_all_comparisons

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

end # end of module