module SigMeasures

export max_percent_increase

"""
    max_percent_increase(results::Dict{Tuple{String, String}, Vector{Float64}}, T::Vector{Int64}, n::Int)

Calculate the t value at which the maximum percent increase over all T values is achieved given the results from all comparisons and return the n key pairs with this t value.

# Arguments
- `results::Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary of results from `k1d_all_comparisons`.
- `T::Vector{Int64}`: A vector of values used in the comparisons.
- `n::Int`: The number of key pairs to return.

# Returns
- `::Vector{Tuple{Tuple{String, String}, Int64}}`: A vector of tuples where each tuple contains a key pair and the t value at which the maximum percent increase is achieved.
"""
function max_percent_increase(results::Dict{Tuple{String, String}, Vector{Float64}}, T::Vector{Int64}, n::Int)::Vector{Tuple{Tuple{String, String}, Int64}}
    theoretical_values = 2 .* T
    max_increases = Vector{Tuple{Tuple{String, String}, Int64}}()

    # Iterate over each key pair in the results dictionary
    for (key_pair, observed_values) in results
        percent_increases = [abs.(observed_values[i] - theoretical_values[i]) / theoretical_values[i] for i in 1:length(T)]
        # Find the index of the maximum percent increase
        max_idx = argmax(percent_increases)
        # Append the key pair and the t value corresponding to the maximum percent increase
        push!(max_increases, (key_pair, T[max_idx]))
    end

    # Sort the max_increases vector by the t value in descending order of the percent increase at that t
    sorted_increases = sort(max_increases, by = x -> x[2], rev = true)

    # Return the top n key pairs with the largest percent increase at their respective t
    return sorted_increases[1:min(n, length(sorted_increases))]
end

end # end of module
