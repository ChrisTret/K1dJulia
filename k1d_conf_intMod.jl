module K1dConfInt

using Statistics

"""
    bootstrap_variance_k1d(Data::Dict{String, Vector{Int64}}, T::Vector{Int64}, n::Int, B::Int) -> Dict{Tuple{String, String}, Vector{Float64}}

Calculate the bootstrap variance of the estimated K function.

# Arguments
- `Data::Dict{String, Vector{Int64}}`: A dictionary where each key is a name and the value is a vector of `genoCenter` values.
- `T::Vector{Int64}`: A vector of values for the diameters.
- `n::Int`: The sample size for each bootstrap sample.
- `B::Int`: The number of bootstrap samples.

# Returns
- `Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary where each key is a tuple of names from `Data` and the value is the bootstrap variance vector.
"""
function bootstrap_variance_k1d(Data::Dict{String, Vector{Int64}}, T::Vector{Int64}, n::Int, B::Int)::Dict{Tuple{String, String}, Vector{Float64}}
    results = Dict{Tuple{String, String}, Vector{Float64}}()
    dict_keys = collect(keys(Data))
    
    # Initialize the parallel processing
    Threads.@threads for key1 in dict_keys
        for key2 in dict_keys
            # Get the genoCenter vectors for the current pair of keys
            X = Data[key1]
            Y = Data[key2]
            K_bootstrap = zeros(Float64, length(T), B)
            
            for b in 1:B
                sample_X = rand(X, n)
                if key1 == key2
                    K_bootstrap[:, b] = k1d_univ(sample_X, T)
                else
                    sample_Y = rand(Y, n)
                    K_bootstrap[:, b] = k1d_biv(sample_X, sample_Y, T)
                end
            end
            
            # Compute the variance for each t in T
            variance_K = var(K_bootstrap, dims=2)
            results[(key1, key2)] = vec(variance_K)
        end
    end
    
    return results
end

end # end of module