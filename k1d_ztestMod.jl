module K1dZTests

using Distributions

export z_test, z_tests_all_comparisons

function z_test(results::Dict{Tuple{String, String}, Vector{Float64}}, T::Vector{Int64}, key_pair::Tuple{String, String}, variance::Vector{Float64}, alpha::Float64)
    μ = 2 .* T
    mean_centered = results[key_pair] .- μ
    z_obs = mean_centered ./ sqrt.(variance)

    z_test = quantile(Normal(), 1 - alpha/2)

    max_deviation_distance = T[argmax(z_obs)]

    clustering = z_obs .> z_test
    dispersion = z_obs .< z_test
    significant = clustering .| dispersion
    
    clustering_distances = T[clustering]
    dispersion_distances = T[dispersion]
    significant_distances = T[significant]
    
    summary = (
        key_pair = key_pair,
        max_deviation_distance = max_deviation_distance,
        clustering_distances = clustering_distances,
        dispersion_distances = dispersion_distances,
        significant_distances = significant_distances,
        z_obs = z_obs
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


function z_tests_all_comparisons(results::Dict{Tuple{String, String}, Vector{Float64}}, T::Vector{Int64}, variances::Dict{Tuple{String, String}, Vector{Float64}}, alpha::Float64)
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

end # end of module