module K1dPlots

export k_plot

using Plots

"""
    k_plot

Plots both the observed K̂ and the theoretical K under CSR against T

# Arguments
`K_dict::Dict{Tuple{String, String}, Vector{Float64}}`: Dictionary containing the data vectors
`keys::Tuple{String, String}`: Tuple of keys to access the data
`T::Vector{Int64}`: Vector of t values searched over

"""
function k_plot(K_dict::Dict{Tuple{String, String}, Vector{Float64}}, keys::Tuple{String, String}, T::Vector{Int64})
    K̂ = K_dict[keys]
    key1, key2 = keys
    
    # E(K̂) under CSR
    theoretical = 2 * T
    
    # Generate title based on keys
    title = "Observed vs Theoretical K Values: $key1 - $key2"
    
    plot(T, theoretical, color = :blue, label = "Theoretical", title = title)
    plot!(T, K̂, color = :red, label = "Observed")
end

end # end of module
