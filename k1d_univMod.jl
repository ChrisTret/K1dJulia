module K1dUniv

export k1d_univ

"""
    k1d_univ(X::Vector{Int64}, T::Vector{Int64})

Estimate Ripley's K function for a univariate point pattern of integer values.

# Arguments
- `X::Vector{Int64}`: Vector with data of one kind (integer values) sorted in ascending order.
- `T::Vector{Int64}`: Vector with diameters (integer values) to search over in ascending order.

# Returns
- `::Vector{Float64}`: Vector of estimated K values for each `t` in `T`.
"""
function k1d_univ(X::Vector{Int64}, T::Vector{Int64})
    N = length(X)
    X_sorted = sort(X)
    
    lambda_hat = N / (X_sorted[end] - X_sorted[1] + 1)
    obs_sum = zeros(Int, length(T))
    K_hat = zeros(Float64, length(T))
    
    for a in 1:N
        down_cum = zeros(Int, length(T))
        up_cum = zeros(Int, length(T))
        down_next_start = 1
        up_next_start = 1
        curr_sum = zeros(Int, length(T))

        for t_idx in 1:length(T)
            t = T[t_idx]
            down_result = check_range_univ(X_sorted, a, down_next_start, t, :down)
            up_result = check_range_univ(X_sorted, a, up_next_start, t, :up)
            
            # Update cumulative sums for the current t
            down_cum[t_idx] = down_result.my_sum
            up_cum[t_idx] = up_result.my_sum
            
            # Update the next start positions
            down_next_start = down_result.next_start
            up_next_start = up_result.next_start
            
            # Add cumulative sum from previous t to current t
            if t_idx > 1
                curr_sum[t_idx] += curr_sum[t_idx - 1]
            end
            
            # Add the current observed values to the obs_sum
            curr_sum[t_idx] += down_cum[t_idx] + up_cum[t_idx]
            obs_sum[t_idx] += curr_sum[t_idx]
        end
    end
    
    K_hat .= obs_sum / (lambda_hat * N)
    return K_hat
end

"""
    check_range_univ(X::Vector{Int64}, a::Int64, start::Int64, t::Int64, direction::Symbol)

Count the number of observations within distance `t` of element `X[a]` in the specified direction.

# Arguments
- `X::Vector{Int64}`: Vector with data of one kind (integer values) sorted in ascending order.
- `a::Int64`: Index of the element in `X` to search around.
- `start::Int64`: Starting index for the search.
- `t::Int64`: Distance to search over.
- `direction::Symbol`: Direction of the search (`:up` or `:down`).

# Returns
- `my_sum::Int64`: Number of observations within `t` of element `X[a]` in the specified direction.
- `next_start::Int64`: Index to start at for the next call.
"""
function check_range_univ(X::Vector{Int64}, a::Int64, start::Int64, t::Int64, direction::Symbol)
    next_start = start
    my_sum = 0
    if direction == :up
        while a + next_start <= length(X) && X[a + next_start] < X[a] + t
            my_sum += 1
            next_start += 1
        end
    elseif direction == :down
        while a > next_start && X[a] - t < X[a - next_start]
            my_sum += 1
            next_start += 1
        end
    end
    return (my_sum = my_sum, next_start = next_start)
end

end # end of Module