module K1dFun

using Base.Threads, Statistics
export k1d_univ, k1d_biv, k1d_all_comparisons, k1d_all_comparisons_indiv_chrom, k1d_mean_across_chromosomes, k1d_mean_across_regions


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
    
    lambda_hat = N / (X[N] - X[1] + 1)
    obs_sum = zeros(Int, length(T))
    K_hat = zeros(Float64, length(T))

    curr_sum = zeros(Int, length(T))
    
    for a in 1:N
        fill!(curr_sum,0)
        down_next_start = 1
        up_next_start = 1

        for t_idx in 1:length(T)
            t = T[t_idx]

            down_sum, down_next_start = check_range_univ(X, a, down_next_start, t, :down)
            up_sum, up_next_start = check_range_univ(X, a, up_next_start, t, :up)
            
            # Add cumulative sum from previous t to current t
            if t_idx > 1
                curr_sum[t_idx] = curr_sum[t_idx - 1]
            end
            
            # Add the current observed values to the obs_sum
            # curr_sum: number around this specific a in X at this specific t
            curr_sum[t_idx] += down_sum + up_sum

            # obs_sum: number around all a in X so far at this specific t
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
   # println("Obs_sum: ", obs_sum)
    K_hat = obs_sum ./ (lambda_hat_x * lambda_hat_y * area)
    
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

    # Loop over each pair of keys (enforce alphabetical order)
    for key1 in dict_keys
        for key2 in dict_keys
            # Skip duplicate pairs by enforcing alphabetical order
            if key1 > key2
                continue
            end

            # Get the genoCenter vectors for the current pair of keys
            X = Data[key1]
            Y = Data[key2]

            if key1 == key2
                result = k1d_univ(X, T)
            else
                result = k1d_biv(X, Y, T)
            end
            # Store results in alphabetical order
            results[(key1, key2)] = result
        end
    end
    return results
end



"""
    k1d_all_comparisons(data::Dict{String, Dict{Int, Vector{Int64}}}, T::Vector{Int64}) 
    -> Dict{Tuple{String, String, Int}, Vector{Float64}}

Computes k1d statistics for all pairs of outer and inner keys in the provided data dictionary, producing either univariate or bivariate results based on the keys.

# Arguments
- `data::Dict{String, Dict{Int, Vector{Int64}}}`: A dictionary where each outer key maps to another dictionary. The inner dictionary maps integer keys (representing regions or categories) to vectors of integer values.
- `T::Vector{Int64}`: A vector of integer values representing the points at which the k1d statistics will be computed.

# Returns
- `Dict{Tuple{String, String, Int}, Vector{Float64}}`: A dictionary where each key is a tuple containing two outer keys and an inner key from the `data` dictionary, and each value is a vector of computed k1d statistics.

# Notes
- For each pair of outer keys (`key1`, `key2`), the function computes either univariate (`k1d_univ`) or bivariate (`k1d_biv`) statistics for matching inner keys. If both vectors are empty, the result for that key pair is a vector of zeros.
- Univariate statistics (`k1d_univ`) are computed if `key1` and `key2` are the same; otherwise, bivariate statistics (`k1d_biv`) are computed.
"""
function k1d_all_comparisons(data::Dict{String, Dict{Int, Vector{Int64}}}, T::Vector{Int64})::Dict{Tuple{String, String, Int}, Vector{Float64}}
    # Initialize the results dictionary
    results = Dict{Tuple{String, String, Int}, Vector{Float64}}()

    # Get the keys of the outer dictionary
    outer_keys = collect(keys(data))

    # Loop over each pair of outer keys
    for key1 in outer_keys
        for key2 in outer_keys
            # Skip duplicate pairs by enforcing alphabetical order
            if key1 > key2
                continue
            end

            # Get the inner dictionaries for the current pair of keys
            inner_dict1 = data[key1]
            inner_dict2 = data[key2]

            # Get the keys of the inner dictionaries
            inner_keys1 = collect(keys(inner_dict1))
            inner_keys2 = collect(keys(inner_dict2))

            # Loop over each pair of matching inner keys
            for inner_key in intersect(inner_keys1, inner_keys2)
                # Get the genoCenter vectors for the current pair of inner keys
                X = inner_dict1[inner_key]
                Y = inner_dict2[inner_key]

                if isempty(X) || isempty(Y)
                    # If either vector is empty, assign a value of 0 for each T
                    results[(key1, key2, inner_key)] = fill(0.0, length(T))
                else
                    if key1 == key2
                        result = k1d_univ(X, T)
                    else
                        result = k1d_biv(X, Y, T)
                    end
                    results[(key1, key2, inner_key)] = result
                end
            end
        end
    end

    return results
end

"""
    k1d_all_comparisons(data::Dict{String, Dict{String, Vector{Int64}}}, T::Vector{Int64}) 
    -> Dict{Tuple{String, String, String}, Vector{Float64}}

Computes k1d statistics for all pairs of outer keys and matching inner keys in the provided data dictionary. The function produces either univariate or bivariate results based on the keys.

# Arguments
- `data::Dict{String, Dict{String, Vector{Int64}}}`: A dictionary where each outer key maps to another dictionary. The inner dictionary maps string keys (representing chromosomes or categories) to vectors of integer values.
- `T::Vector{Int64}`: A vector of integer values representing the points at which the k1d statistics will be computed.

# Returns
- `Dict{Tuple{String, String, String}, Vector{Float64}}`: A dictionary where each key is a tuple containing two outer keys and an inner key (chromosome) from the `data` dictionary, and each value is a vector of computed k1d statistics.

# Notes
- For each pair of outer keys (`key1`, `key2`), the function computes either univariate (`k1d_univ`) or bivariate (`k1d_biv`) statistics for matching inner keys. 
- Univariate statistics (`k1d_univ`) are computed if `key1` and `key2` are the same; otherwise, bivariate statistics (`k1d_biv`) are computed.
"""
function k1d_all_comparisons(data::Dict{String, Dict{String, Vector{Int64}}}, T::Vector{Int64})::Dict{Tuple{String, String, String}, Vector{Float64}}
    results = Dict{Tuple{String, String, String}, Vector{Float64}}()
    outer_keys = keys(data)

    for key1 in outer_keys
        dict1 = data[key1]
        inner_keys1 = keys(dict1)
        
        for key2 in outer_keys
            # Skip duplicate pairs by enforcing alphabetical order
            if key1 > key2
                continue
            end

            dict2 = data[key2]

            for chrom in inner_keys1
                if haskey(dict2, chrom)
                    X = dict1[chrom]
                    Y = dict2[chrom]

                    if key1 == key2
                        result = k1d_univ(X, T)
                    else
                        result = k1d_biv(X, Y, T)
                    end
                    # Store results in alphabetical order
                    results[(key1, key2, chrom)] = result
                end
            end
        end
    end
    return results
end



"""
    k1d_all_comparisons(data::Dict{String, Dict{String, Dict{Int, Vector{Int64}}}}, T::Vector{Int64}) 
    -> Dict{Tuple{String, String, String, Int}, Vector{Float64}}

Computes k1d statistics for all pairs of outer keys, matching inner keys, and region identifiers in the provided data dictionary. The function produces either univariate or bivariate results based on the keys.

# Arguments
- `data::Dict{String, Dict{String, Dict{Int, Vector{Int64}}}}`: A nested dictionary where each outer key maps to another dictionary, which further maps string keys (representing chromosomes) to yet another dictionary. The innermost dictionary maps integer region identifiers to vectors of integer values.
- `T::Vector{Int64}`: A vector of integer values representing the points at which the k1d statistics will be computed.

# Returns
- `Dict{Tuple{String, String, String, Int}, Vector{Float64}}`: A dictionary where each key is a tuple containing two outer keys, a chromosome, and a region identifier, and each value is a vector of computed k1d statistics.
"""
function k1d_all_comparisons(data::Dict{String, Dict{String, Dict{Int, Vector{Int64}}}}, T::Vector{Int64})::Dict{Tuple{String, String, String, Int}, Vector{Float64}}
    results = Dict{Tuple{String, String, String, Int}, Vector{Float64}}()
    outer_keys = keys(data)

    for key1 in outer_keys
        dict1 = data[key1]
        inner_keys1 = keys(dict1)
        
        for key2 in outer_keys
            # Enforce alphabetical ordering by swapping keys if necessary
            if key1 > key2
                continue
            end

            dict2 = data[key2]

            for chrom in inner_keys1
                if haskey(dict2, chrom)
                    regions1 = keys(dict1[chrom])
                    regions2 = keys(dict2[chrom])

                    for region in regions1
                        if region in regions2
                            X = dict1[chrom][region]
                            Y = dict2[chrom][region]

                            # Check for empty vectors and assign k values of 0 if empty
                            if isempty(X) || isempty(Y)
                                result = fill(0.0, length(T))
                            else
                                if key1 == key2
                                    result = k1d_univ(X, T)
                                else
                                    result = k1d_biv(X, Y, T)
                                end
                            end
                            results[(key1, key2, chrom, region)] = result
                        end
                    end
                end
            end
        end
    end
    return results
end


"""
    k1d_mean_across_chromosomes(k1d_results::Dict{Tuple{String, String, String}, Vector{Float64}})
    -> Dict{Tuple{String, String}, Vector{Float64}}

Calculates the mean of k1d statistics across different chromosomes for each pair of element names.

# Arguments
- `k1d_results::Dict{Tuple{String, String, String}, Vector{Float64}}`: A dictionary where each key is a tuple containing two element names and a chromosome identifier, and each value is a vector of k1d statistics corresponding to that key.

# Returns
- `Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary where each key is a tuple containing two element names, and each value is a vector representing the mean k1d statistics across all chromosomes for that pair of element names.
"""
function k1d_mean_across_chromosomes(k1d_results::Dict{Tuple{String, String, String}, Vector{Float64}})::Dict{Tuple{String, String}, Vector{Float64}}
    # Initialize the dictionary to store mean values
    mean_results = Dict{Tuple{String, String}, Vector{Float64}}()
    
    # Extract unique pairs of element names from k1d_results keys
    pairs = unique((key[1], key[2]) for key in keys(k1d_results))
    
    # Loop over each pair of element names
    for (key1, key2) in pairs
        # Initialize a matrix to store k1d values for all matching chromosomes
        k1d_matrix = []

        # Extract k1d values for all matching chromosomes
        for (key, value) in k1d_results
            if (key[1], key[2]) == (key1, key2)
                push!(k1d_matrix, value)
            end
        end

        if !isempty(k1d_matrix)
            k1d_matrix = reduce(hcat, k1d_matrix)
            
            # Calculate the mean for each t
            mean_values = mean(k1d_matrix, dims=2)
            mean_results[(key1, key2)] = vec(mean_values)
        end
    end
    
    return mean_results
end



"""
    k1d_mean_across_regions(k1d_results::Dict{Tuple{String, String, String, Int}, Vector{Float64}})
    -> Dict{Tuple{String, String}, Vector{Float64}}

Calculates the mean of k1d statistics across different chromosomes and regions for each pair of element names.

# Arguments
- `k1d_results::Dict{Tuple{String, String, String, Int}, Vector{Float64}}`: A dictionary where each key is a tuple containing two element names, a chromosome identifier, and a region identifier, and each value is a vector of k1d statistics corresponding to that key.

# Returns
- `Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary where each key is a tuple containing two element names, and each value is a vector representing the mean k1d statistics across all chromosomes and regions for that pair of element names.
"""
function k1d_mean_across_regions(k1d_results::Dict{Tuple{String, String, String, Int}, Vector{Float64}})::Dict{Tuple{String, String}, Vector{Float64}}
    # Initialize the dictionary to store mean values
    mean_results = Dict{Tuple{String, String}, Vector{Float64}}()
    
    # Group data by element pairs to avoid repetitive filtering
    grouped_results = Dict{Tuple{String, String}, Vector{Vector{Float64}}}()
    for (key, value) in k1d_results
        pair = (key[1], key[2])
        if !haskey(grouped_results, pair)
            grouped_results[pair] = []
        end
        push!(grouped_results[pair], value)
    end

    # Compute mean values for each pair
    for (pair, matrices) in grouped_results
        # Combine all vectors into a single matrix
        k1d_matrix = reduce(hcat, matrices)
        # Calculate the mean for each t
        mean_values = mean(k1d_matrix, dims=2)
        mean_results[pair] = vec(mean_values)
    end
    
    return mean_results
end



"""
    k1d_mean_across_regions(k1d_results::Dict{Tuple{String, String, Int}, Vector{Float64}})
    -> Dict{Tuple{String, String}, Vector{Float64}}

Calculates the mean of k1d statistics across different regions for each pair of element names.

# Arguments
- `k1d_results::Dict{Tuple{String, String, Int}, Vector{Float64}}`: A dictionary where each key is a tuple containing two element names and a region identifier, and each value is a vector of k1d statistics corresponding to that key.

# Returns
- `Dict{Tuple{String, String}, Vector{Float64}}`: A dictionary where each key is a tuple containing two element names, and each value is a vector representing the mean k1d statistics across all regions for that pair of element names.
"""
function k1d_mean_across_regions(k1d_results::Dict{Tuple{String, String, Int}, Vector{Float64}})::Dict{Tuple{String, String}, Vector{Float64}}
    # Initialize the dictionary to store mean values
    mean_results = Dict{Tuple{String, String}, Vector{Float64}}()

    # Extract unique pairs of element names from k1d_results keys
    pairs = unique((key[1], key[2]) for key in keys(k1d_results))
    
    # Loop over each pair of element names
    for (key1, key2) in pairs
        # Initialize a matrix to store k1d values for all matching regions
        k1d_matrix = []

        # Extract k1d values for all matching regions
        for (key, value) in k1d_results
            if (key[1], key[2]) == (key1, key2)
                push!(k1d_matrix, value)
            end
        end

        if !isempty(k1d_matrix)
            k1d_matrix = reduce(hcat, k1d_matrix)
            
            # Calculate the mean for each t
            mean_values = mean(k1d_matrix, dims=2)
            mean_results[(key1, key2)] = vec(mean_values)
        end
    end
    
    return mean_results
end

end # end of module