module K1dPlots

export k_plot, plot_top_n

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


"""
    plot_top_n

Utilizes k_plot to generate individual plots for the top n key pairs, each highlighting the t value where the maximum percent increase occurs.

# Arguments
- `results::Dict{Tuple{String, String}, Vector{Float64}}`: Dictionary containing the data vectors.
- `T::Vector{Int64}`: Vector of t values searched over.
- `top_n::Vector{Tuple{Tuple{String, String}, Int64}}`: Vector of tuples where each tuple contains a key pair and the t value where the maximum percent increase is achieved.
"""
function plot_top_n(results::Dict{Tuple{String, String}, Vector{Float64}}, T::Vector{Int64}, top_n::Vector{Tuple{Tuple{String, String}, Int64}})
    for (keys, max_t) in top_n
        k_plot(results, keys, T)  # Use k_plot to handle individual plotting
        
        # Additional logic to highlight the t value with maximum percent increase (optional)
        K̂ = results[keys]
        max_t_idx = findfirst(==(max_t), T)
        scatter!([T[max_t_idx]], [K̂[max_t_idx]], label="Max Inc at t=$max_t", color=:green, markershape=:circle)
        
        # Customize and display the plot
        title!("$keys, Max Inc at t=$max_t")
        display(current())  # Explicitly display the current plot
    end
end

end # end of module