module KPlots

export k_plot

"""
    k_plot

Plots both the observed K̂ and the theoretical K under CSR against T

# Arguments
`K̂::Vector{Float64}`: Vector of observed K̂ values from k1d_univ
`T::Vector{Int64}`: Vector of t values searched over

"""
function k_plot(K̂::Vector{Float64}, T::Vector{Int64})
   
    # E(K̂) under CSR
    theoretical = 2 * T

    plot(T, theoretical, color = :blue, label = "Theoretical")
    plot!(T, K̂, color = :red, label = "Observed")

end

end # end of module

