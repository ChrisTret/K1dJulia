module K1dFun

export k1d_univ, k1d_biv, k1d_all_comparisons

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

"""
    k1d_biv(X::Vector{Int64}, Y::Vector{Int64}, T::Vector{Int64})

Estimate Ripley's K function for a bivariate point pattern of integer values.

# Arguments
- `X::Vector{Int64}`: Vector with data of one kind (integer values) sorted in ascending order.
- `Y::Vector{Int64}`: Vector with data of another kind (integer values) sorted in ascending order.
- `T::Vector{Int64}`: Vector with diameters (integer values) to search over in ascending order.

# Returns
- `::Vector{Float64}`: Vector of estimated K_XY values for each `t` in `T`.
"""
function k1d_biv(X::Vector{Int64}, Y::Vector{Int64}, T::Vector{Int64})
    Nx = length(X)
    Ny = length(Y)
    
    area = maximum([maximum(X), maximum(Y)]) - minimum([minimum(X), minimum(Y)]) + 1    
    lambda_hat_x = Nx / area
    lambda_hat_y = Ny / area
    
    obs_sum = zeros(length(T))
    K_hat = zeros(length(T))
    
    for t in 1:length(T)
        start = 1
        for a in 1:Nx
            result = check_range_biv(X, Y, a, start, T[t])
            obs_sum[t] += result[1]
            start = result[2]
        end
    end
    
    K_hat = obs_sum ./ (lambda_hat_y * lambda_hat_y * area)
    
    return K_hat
end

"""
    check_range_biv(X::Vector{Int64}, Y::Vector{Int64}, a::Int64, j::Int64, t::Int64)

Count the number of observations in `Y` within distance `t` of element `X[a]`.

# Arguments
- `X::Vector{Int64}`: Vector with data of one kind (integer values) sorted in ascending order.
- `Y::Vector{Int64}`: Vector with data of another kind (integer values) sorted in ascending order.
- `a::Int64`: Index of the element in `X` to search around.
- `j::Int64`: Starting index for the search in `Y`.
- `t::Int64`: Distance to search over.

# Returns
- `my_sum::Int64`: Number of observations in `Y` within `t` of element `X[a]`.
- `next_start::Int64`: Index to start at for the next call in `Y`.
"""
function check_range_biv(X::Vector{Int64}, Y::Vector{Int64}, a::Int64, j::Int64, t::Int64)
    next_start = j
    my_sum = 0
    
    while j <= length(Y) && X[a] + t > Y[j]
        if X[a] - t < Y[j]
            my_sum += 1
        else
            next_start = j + 1
        end
        j += 1
    end
    
    return my_sum, next_start
end


"""
    k1d_biv_all_comparisons(Data::Dict{String, Vector{Int64}}, T::Vector{Int64}) -> Dict{Tuple{String, String}, Vector{Int64}}

Runs the `k1d_biv` function for each possible pair of keys in the dictionary `Data`. The result for each pair is stored in a nested dictionary where the key is a tuple of the form `(key1, key2)` and the value is the result vector from `k1d_biv`.

# Arguments
- `Data::Dict{String, Vector{Int64}}`: A dictionary where each key is a name (e.g., `repFamily`) and the value is a vector of `genoCenter` values.
- `T::Vector{Int64}`: A vector of values to be used as the third parameter for the `k1d_biv` function.

# Returns
- `Dict{Tuple{String, String}, Vector{Int64}}`: A dictionary where each key is a tuple of two names from `Data` and the value is the result vector from `k1d_biv` for the corresponding pair of `genoCenter` vectors.

"""
function k1d_all_comparisons(Data::Dict{String, Vector{Int64}}, T::Vector{Int64})::Dict{Tuple{String, String}, Vector{Float64}}
    # Initialize the results dictionary
    results = Dict{Tuple{String, String}, Vector{Float64}}()
    # Get the keys of the dictionary
    dict_keys = collect(keys(Data))


# Loop over each pair of keys (including the same-key pair)
    for key1 in dict_keys
        for key2 in dict_keys
            # Get the genoCenter vectors for the current pair of keys
            X = Data[key1]
            Y = Data[key2]

            if key1 == key2
                result = k1d_univ(X, T)
            else
                result = k1d_biv(X, Y, T)
            end
            results[(key1, key2)] = result
        end
    end
return results
end

end # end of module
