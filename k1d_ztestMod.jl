module K1dZTests

using Distributions, Printf

export z_test, z_tests_all_comparisons, print_z_summary, z_proportions, top_n_means


function z_to_p(z_score::Float64)::Float64
    dist = Normal(0, 1)
    p_value = 2 * (1 - cdf(dist, abs(z_score))) # Two-tailed p-value
    return p_value
end


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

    if any(significant)
        println(key_pair," ", max_deviation_distance)
    end

    # println("Clustering detected at t = ", T[clustering])
    # println("Dispersion detected at t = ", T[dispersion])
    # println("Significant deviation from CSR detected at t = ", T[significant])
    # println("Most significant distance ", max_deviation_distance)

end


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

    # println("\nSignificant Distances:")
    # if !isempty(summary.significant_distances)
    #     println(join([@sprintf("%5d", d) for d in summary.significant_distances], " "))
    # else
    #     println("None")
    # end

    println("\nZ-Statistics:")
    println(join([format_z_score(z) for z in summary.z_scores], " "))
    
    println("\n")
end


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