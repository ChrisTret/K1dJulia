
# Demo with length 10 vec, t=3
#demo = 1:5:100
#k1d_univ(demo, T = [10])

my_test = [1, 6, 12, 15, 17, 24]

# Algorithm idea
# 1. select point a
# 2. calculate values in (a-t, a)
# 3. calculated values in (a, a+t)
# 4. sum values
# 5. repeat for next i

################################################################################
# Estimate Ripley's K function with loops
#
# @param X vector with data of one kind sorted in ascending order
# @param T vector with diameters to search over in ascending order
#
# @returns vector of estimated K value for each t in T
################################################################################
function k1d_univ(X, T)
    N = length(X)
    X = sort(X)
    
    # println("Obs X: $N")
    # println("Distance: $(X[end] - X[1] + 1)")
    
    lambda_hat = N / (X[end] - X[1] + 1)  # number occurrences over length of interval
    # println("lambda: $lambda_hat")
    
    down_sum = zeros(Int, length(T))
    up_sum = zeros(Int, length(T))
    obs_sum = zeros(Int, length(T))
    K_hat = zeros(Float64, length(T))
    
    for a in 1:N
        for t in 1:length(T)
            down_sum[t] = check_interval_down(X, a, T[t])
            up_sum[t] = check_interval_up(X, a, T[t])
            obs_sum[t] += down_sum[t] + up_sum[t]
        end
    end

    K_hat .= obs_sum / (lambda_hat * N)

    return K_hat
end

################################################################################
# Count number of observations within t of element a in X on the right
#
# @param X vector with data of one kind sorted in ascending order
# @param a index of element in X to search around
# @param t length to search over
#
# @returns number of observations within t of element a in X on the right
################################################################################
function check_interval_up(X, a, t)
    j = 1
    my_sum = 0
    while a + j <= length(X) && X[a + j] < X[a] + t
        my_sum += 1
        j += 1
    end
    return my_sum
end

################################################################################
# Count number of observations within t of element a in X on the left
#
# @param X vector with data of one kind sorted in ascending order
# @param a index of element in X to search around
# @param t length to search over
# 
# @returns number of observations within t of element a in X on the left
################################################################################
function check_interval_down(X, a, t)
    j = 1
    my_sum = 0
    while a > j && X[a] - t < X[a - j]
        my_sum += 1
        j += 1
    end
    return my_sum
end

################################################################################
# Compare speed of two Ripley's K algorithms
#
# @param X vector with data of one kind sorted in ascending order
# @param T vector with diameters to search over in ascending order
# @param function1 first function to compare
# @param function2 second function to compare
# @param repetitions number of times to run each experiment (default 50)
#
# @returns vector of estimated K value for each t in T
################################################################################
using Statistics, Plots

function speed_comparison(X, T, function1, function2, repetitions=50)
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
    
    hist1 = histogram(times1, alpha=0.6, label="Function 1", legend=:topright)
    histogram!(hist1, times2, alpha=0.6, label="Function 2")
    
    display(hist1)
    
    return ttest_ind(times1, times2)
end

################################################################################
# Calculate number of data points within t of edge
#
# @param X vector with data of one kind sorted in ascending order
# @param t diameter to check for data within t of edge
#
# @returns number of data points within t of edge
################################################################################
function danger_zone(X, t)
    danger_count1 = 1  # X[1] is by definition in the danger zone
    danger_count2 = 1  # X[end] is by definition in the danger zone
    n = length(X)
    
    while X[danger_count1 + 1] - t < X[1]
        danger_count1 += 1
    end
    while X[n - danger_count2] + t > X[end] && n - danger_count2 > danger_count1 + 1
        danger_count2 += 1
    end
    
    danger_count = danger_count1 + danger_count2
    
    return danger_count
end

# # Example usage
# demo = collect(10:1000)
# k1d_univ(demo, [10])

my_test = [1, 6, 12, 15, 17, 24]
k1d_univ_optimized(my_test, [3, 4, 10])


#elapsed = @elapsed k1d_univ(demo, [3, 4, 10]) # hand-calculated 4/3, 8/3. 32/3

# k=5
# elapsed = zeros(Float64,k)
# for i in 1:k
#     top = 10 * 10^i
#     demo = 1:top
#     elapsed[i] = @elapsed k1d_univ(demo, 1000:500:5000) # hand-calculated 4/3, 8/3. 32/3
# end
# println(elapsed)
