module ZTests

using Distributions

export z_test

function z_test(results::Dict{Tuple{String, String}, Vector{Float64}}, T::Vector{Int64}, key_pair::Tuple{String, String}, variance::Vector{Float64}, alpha::Float64)
    μ = 2 * T
    mean_centered = results[key_pair] .- μ
    z_obs = mean_centered ./ sqrt.(variance)

    z_test = quantile(Normal(), 1 - alpha/2)

    max_deviation_distance = T[argmax(z_obs)]

    clustering = z_obs .> z_test
    dispersion = z_obs .< z_test
    significant = clustering .| dispersion
    

    println("Clustering detected at t = ", T[clustering])
    println("Dispersion detected at t = ", T[dispersion])
    println("Significant deviation from CSR detected at t = ", T[significant])
    println("Most significant distance", max_deviation_distance)

end

end # end of module