module SpeedTests

export speed_comparison_univ

#check if this is the rigth package
using Plots 

# Compare speed of the original and optimized functions
function speed_comparison_univ(X::Vector{Int64}, T::Vector{Int64}, function1::Function, function2::Function, repetitions=50::Int64, title = "Speed Comparison")
    println("Starting Speed Comparison")
    total_time1 = zeros(repetitions)
    total_time2 = zeros(repetitions)
    
    for i in 1:repetitions
        start = time()
        function1(X, T)
        total_time1[i] = time() - start
    end
    
    for i in 1:repetitions
        start = time()
        function2(X, T)
        total_time2[i] = time() - start
    end
    
    times1 = total_time1
    times2 = total_time2
    
    combined_times = vcat(times1, times2)
    min_time = minimum(combined_times)
    max_time = maximum(combined_times)
    bin_edges = range(min_time, 0.004, step=0.0005)


    hist1 = histogram(times1, bins = 40, alpha=0.6, label="Function 1", legend=:topright, title = title)
    histogram!(hist1, times2, bins = 40, alpha=0.6, label="Function 2")
    
    display("image/png", hist1)

    println("Mean time for Function 1: ", mean(times1))
    println("Mean time for Function 2: ", mean(times2))
    println("Total Simulation Time: ", sum(total_time1) + sum(total_time2))


    return nothing
end

end # end of module