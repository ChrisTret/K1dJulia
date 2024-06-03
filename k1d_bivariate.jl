module K1dBiv

export k1d_biv

"""
    ripley_k(X, Y, T)

Estimate Ripley's K function for a bivariate point pattern.

# Arguments
- `X::AbstractVector`: Vector with data of one kind sorted in ascending order.
- `Y::AbstractVector`: Vector with data of one kind sorted in ascending order.
- `T::AbstractVector`: Vector with diameters to search over in ascending order.

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

end # end of module