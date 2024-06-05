module KPlots

export K_plot
function K_plot(K̂::Vector{Float64}, T::Vector{Int64})
   
    # E(K̂) under CSR
    theoretical = 2 * T

    plot(T, theoretical, color = :blue, label = "Theoretical")
    plot!(T, K̂, color = :red, label = "Observed")

end


end

