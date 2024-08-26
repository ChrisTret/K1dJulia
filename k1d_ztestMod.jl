module K1dZTests

using Distributions, Printf

export z_test, z_tests_all_comparisons, print_z_summary, z_proportions, top_n_means


"""
    z_to_p(z_score::Float64) -> Float64

Calculates the two-tailed p-value for a given z-score, assuming a standard normal distribution.

# Arguments
- `z_score::Float64`: The z-score for which to calculate the p-value.

# Returns
- `Float64`: The two-tailed p-value corresponding to the input z-score.
"""
function z_to_p(z_score::Float64)::Float64
    dist = Normal(0, 1)
    p_value = 2 * (1 - cdf(dist, abs(z_score))) # Two-tailed p-value
    return p_value
end

"""
    z_test(results::Dict{Tuple{String, String}, Vector{Float64}}, T::Vector{Int64}, 
           key_pair::Tuple{String, String}, variance::Vector{Float64}, 
           alpha::Float64 = 0.05) -> NamedTuple

Performs a z-test on the provided data to assess significant deviations from the expected mean based on the given variances.

# Arguments
- `results::Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary containing vectors of results for each key pair.
- `T::Vector{Int64}`: A vector of distances or time points at which tests are performed.
- `key_pair::Tuple{String, String}`: A tuple representing the key pair for which the z-test is performed.
- `variance::Vector{Float64}`: A vector of variances corresponding to the result vectors for the given key pair.
- `alpha::Float64 = 0.05`: Significance level for the test (default is 0.05).

# Returns
- `NamedTuple`: A summary of the z-test results, including z-scores, the most significant distance and its z-score, top distances with their corresponding z-scores and p-values, clustering and dispersion distances, significant distances, and the proportion of z-scores above the significance threshold.
"""
function z_test(results::Dict{Tuple{String, String}, Vector{Float64}}, T::Vector{Int64}, key_pair::Tuple{String, String}, variance::Vector{Float64}, alpha::Float64 = 0.05)
    μ = 2 .* T
    mean_centered = results[key_pair] .- μ
    z_obs = mean_centered ./ sqrt.(variance)
    p_vals = z_to_p.(z_obs)

    z_test = quantile(Normal(), 1 - alpha/2)
    max_abs_z = maximum(abs.(z_obs))
    
    if max_abs_z > z_test
        most_sig_distance = T[argmax(abs.(z_obs))]
        most_sig_distance_z = z_obs[argmax(abs.(z_obs))] # retains Z-score sign

        # Find the largest z-scores in absolute value and filter based on z_test
        sorted_indices = sortperm(abs.(z_obs), rev=true)
        top_indices = filter(i -> abs(z_obs[i]) > z_test, sorted_indices)
        top_indices = first(top_indices, 10)
        top_distances = T[top_indices]
        top_z_scores = z_obs[top_indices]
        top_p_vals = p_vals[top_indices]
    else
        most_sig_distance = nothing
        most_sig_distance_z = nothing
        top_distances = nothing
        top_z_scores = nothing
        top_p_vals = nothing
    end

    clustering = z_obs .> z_test
    dispersion = z_obs .< -z_test
    significant = clustering .| dispersion
    
    clustering_distances = T[clustering]
    dispersion_distances = T[dispersion]
    significant_distances = T[significant]

    proportion_above_threshold = mean(abs.(z_obs) .> z_test)

    summary = (
        key_pair = key_pair,
        z_scores = z_obs,
        most_sig_distance = most_sig_distance,
        most_sig_distance_z = most_sig_distance_z,
        top_distances = top_distances,
        top_z_scores = top_z_scores,
        top_p_vals = top_p_vals,
        clustering_distances = clustering_distances,
        dispersion_distances = dispersion_distances,
        significant_distances = significant_distances,
        proportion_above_threshold = proportion_above_threshold
    )
    
    return summary
end



"""
    z_tests_all_comparisons(results::Dict{Tuple{String, String}, Vector{Float64}}, 
                            T::Vector{Int64}, variances::Dict{Tuple{String, String}, Vector{Float64}}, 
                            alpha::Float64 = 0.05) -> Dict{Tuple{String, String}, NamedTuple}

Performs z-tests for all key pairs in the `results` and `variances` dictionaries and returns a summary of the test results for each key pair.

# Arguments
- `results::Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary where each key pair maps to a vector of result values.
- `T::Vector{Int64}`: A vector of distances or time points at which tests are performed.
- `variances::Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary where each key pair maps to a vector of variances corresponding to the result vectors in `results`.
- `alpha::Float64 = 0.05`: The significance level for the z-tests (default is 0.05).

# Returns
- `Dict{Tuple{String, String}, NamedTuple}`: A dictionary where each key pair is associated with a summary of z-test results. The summary includes z-scores, the most significant distance and its z-score, top distances with their corresponding z-scores and p-values, clustering and dispersion distances, significant distances, and the proportion of z-scores above the significance threshold.

# Notes
- If a key pair in `results` does not have a corresponding variance in `variances`, the function prints a message indicating the missing variance and skips the z-test for that key pair.
"""
function z_tests_all_comparisons(results::Dict{Tuple{String, String}, Vector{Float64}}, T::Vector{Int64}, variances::Dict{Tuple{String, String}, Vector{Float64}}, alpha::Float64 = 0.05)
    summaries = Dict{Tuple{String, String}, NamedTuple}()

    for key_pair in keys(results)
        if haskey(variances, key_pair)
            summaries[key_pair] = z_test(results, T, key_pair, variances[key_pair], alpha)
        else
            println("Variance not found for key pair: ", key_pair)
        end
    end

    return summaries
end



"""
    format_z_score(z::Float64) -> String

Formats a z-score for display by capping extreme values at -10 and 10 and formatting within a specific range.

# Arguments
- `z::Float64`: The z-score to format.

# Returns
- `String`: A formatted string representation of the z-score. 
    - Returns `">10"` if the z-score is greater than 10.
    - Returns `"<-10"` if the z-score is less than -10.
    - Returns a string with two decimal places for z-scores within the range [-10, 10].
"""
function format_z_score(z::Float64)::String
    z = min(max(z, -10.0), 10.0) # Cap the z-score between -10 and 10
    if z == 10.0
        return ">10"
    elseif z == -10.0
        return "<-10"
    else
        return @sprintf("% .2f", z)
    end
end



"""
    print_z_summary(summary::NamedTuple)

Prints a formatted summary of z-test results contained in the provided `summary`.

# Arguments
- `summary::NamedTuple`: A named tuple containing the results of a z-test, including information such as the key pair, z-scores, most significant distance, top distances with their corresponding z-scores and p-values, clustering distances, and dispersion distances.

# Returns
- `Nothing`: This function prints the summary directly to the console and does not return any value.
"""
function print_z_summary(summary::NamedTuple)
    println("\nSummary for key pair: ", summary.key_pair)

    if summary.most_sig_distance !== nothing
        println("\nMost Significant Distance: ", summary.most_sig_distance)
        println("Associated Z-score: ", summary.most_sig_distance_z)
        if abs(summary.most_sig_distance_z) > 10
            println(
                "Warning: Unreasonably large Z-scores are likely indicators that the associated distance " *
                "is not long enough to reach sufficiently outside the element (gene, TE, TF, etc.). " *
                "Please see the Top Z-scores along with other results below to get a more holist view of the situation."
            )
        end
    else
        println("\nMost Significant Distance: None")
        println("Associated Z-score: None")
    end

    if summary.top_z_scores !== nothing && summary.top_p_vals !== nothing
        println("\nTop Z-scores:")
        println(@sprintf("%-10s %-10s %-10s", "Distance", "Z-score", "P-value"))
        for (distance, z_score, p_value) in zip(summary.top_distances, summary.top_z_scores, summary.top_p_vals)
            z_str = format_z_score(z_score)
            println(@sprintf("%5d %10s %10.4f", distance, z_str, p_value))
        end
    end
    println("\nClustering Distances:")
    if summary.clustering_distances !== nothing
        println(join([@sprintf("%5d", d) for d in summary.clustering_distances], " "))
    else
        println("None")
    end

    println("\nDispersion Distances:")
    if summary.dispersion_distances !== nothing
        println(join([@sprintf("%5d", d) for d in summary.dispersion_distances], " "))
    else
        println("None")
    end

    println("\nZ-Statistics:")
    println(join([format_z_score(z) for z in summary.z_scores], " "))
    
    println("\n")
end



"""
    z_proportions(summaries::Dict{Tuple{String, String}, NamedTuple}, n::Int=10) -> Vector{Tuple{Tuple{String, String}, Float64}}

Calculates and prints the top `n` key pairs with the largest proportions of significant z-scores above a given threshold from the provided summaries.

# Arguments
- `summaries::Dict{Tuple{String, String}, NamedTuple}`: A dictionary where each key pair is associated with a summary of z-test results, including the proportion of z-scores above the significance threshold.
- `n::Int=10`: The number of top key pairs to display based on their proportions of significant z-scores (default is 10).

# Returns
- `Vector{Tuple{Tuple{String, String}, Float64}}`: A vector containing the top `n` key pairs and their corresponding proportions of significant z-scores above the threshold, sorted in descending order.
"""
function z_proportions(summaries::Dict{Tuple{String, String}, NamedTuple}, n::Int=10)
    # Create a vector to hold (key, proportion) tuples
    proportions = [(key, summaries[key].proportion_above_threshold) for key in keys(summaries)]
    
    # Sort the vector by proportions in descending order
    sorted_proportions = sort(proportions, by = x -> x[2], rev = true)

    # Get the top n key pairs and their proportions
    top_n_keys_and_proportions = first(sorted_proportions, n)

    # Print the result in the desired format
    println("Pairs with largest significance proportions")
    println(@sprintf("%-15s\t%-50s", "Sig. Proportion", "Pairs"))
    println("=" ^ 65)
    for (pair, proportion) in top_n_keys_and_proportions
        println(@sprintf("%-15.4f\t%-50s", proportion, pair))
    end

    return top_n_keys_and_proportions
end



"""
    top_n_means(z_scores_dict::Dict{Tuple{String, String}, NamedTuple}, n::Int=10) 
    -> Tuple{Vector{Tuple{Tuple{String, String}, Float64}}, Vector{Tuple{Tuple{String, String}, Float64}}}

Calculates and prints the top `n` and bottom `n` key pairs with the largest and smallest mean z-scores, respectively, from the provided z-scores dictionary.

# Arguments
- `z_scores_dict::Dict{Tuple{String, String}, NamedTuple}`: A dictionary where each key pair is associated with a summary of z-test results, including a vector of z-scores.
- `n::Int=10`: The number of top and bottom key pairs to display based on their mean z-scores (default is 10).

# Returns
- `Tuple{Vector{Tuple{Tuple{String, String}, Float64}}, Vector{Tuple{Tuple{String, String}, Float64}}}`: 
  A tuple containing two vectors:
  1. The first vector contains the top `n` key pairs and their corresponding mean z-scores, sorted in descending order.
  2. The second vector contains the bottom `n` key pairs and their corresponding mean z-scores, sorted in ascending order.
"""
function top_n_means(z_scores_dict::Dict{Tuple{String, String}, NamedTuple}, n::Int=10)
    # Create a vector to hold (key, mean_z) tuples
    mean_z_scores = [(key, mean(value.z_scores)) for (key, value) in z_scores_dict]
    
    # Sort the vector by mean Z-scores in descending order
    sorted_means = sort(mean_z_scores, by = x -> x[2], rev = true)

    # Get the top n key pairs and their mean Z-scores
    top_n_keys_and_means = first(sorted_means, n)

    # Get the bottom n key pairs and their mean Z-scores
    bottom_n_keys_and_means = last(sorted_means, n)

    # Print the result in the desired format
    println("Pairs with largest mean Z-scores")
    println(@sprintf("%-15s\t%-50s", "Mean Z-score", "Pairs"))
    println("=" ^ 65)
    for (pair, mean_z) in top_n_keys_and_means
        println(@sprintf("%-15.4f\t%-50s", mean_z, pair))
    end

    println("\nPairs with smallest mean Z-scores")
    println(@sprintf("%-15s\t%-50s", "Mean Z-score", "Pairs"))
    println("=" ^ 65)
    for (pair, mean_z) in bottom_n_keys_and_means
        println(@sprintf("%-15.4f\t%-50s", mean_z, pair))
    end

    return (top_n_keys_and_means, bottom_n_keys_and_means)
end


end # end of module