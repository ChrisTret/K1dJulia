module K1dHeatmap

using Plots, Clustering, CSV, DataFrames, Distances, NamedArrays

export zscore_heatmap, extract_most_sig_z, create_matrix, export_matrix_to_csv, transform_to_dist_matrix, plot_heatmap_with_clustering, create_effect_size_csv


function extract_most_sig_z(summaries::Dict{Tuple{String, String}, NamedTuple})
    sig_z_dict = Dict{Tuple{String, String}, Float64}()
    for (key_pair, summary) in summaries
        if !isnothing(summary.most_sig_distance_z)
            sig_z_dict[key_pair] = summary.most_sig_distance_z
        else
            sig_z_dict[key_pair] = NaN
        end
    end
    return sig_z_dict
end


function extract_unique_keys(sig_z_dict::Dict{Tuple{String, String}, Float64})
    # Extract the unique keys from the dictionary keys
    keys_flattened = vcat([key[1] for key in keys(sig_z_dict)], [key[2] for key in keys(sig_z_dict)])
    unique_keys = unique(keys_flattened)  # Remove duplicates and return unique keys
    return unique_keys
end


function create_matrix(sig_z_dict::Dict{Tuple{String, String}, Float64})
    # Extract unique keys from the dictionary
    unique_keys = extract_unique_keys(sig_z_dict)
    
    # Sort unique keys alphabetically or by any preferred criteria
    sorted_keys = sort(collect(unique_keys))
    println(sorted_keys)
    
    n = length(sorted_keys)
    
    # Initialize the heatmap matrix with NaN values
    heatmap_matrix = fill(NaN, n, n)

    # Fill the matrix with values from sig_z_dict, matching the sorted order of keys
    for ((key1, key2), z_val) in sig_z_dict
        i = findfirst(isequal(key1), sorted_keys)  # Find index of key1 in sorted_keys
        j = findfirst(isequal(key2), sorted_keys)  # Find index of key2 in sorted_keys
        if !isnothing(i) && !isnothing(j)
            heatmap_matrix[i, j] = z_val
        end
    end
    println(typeof(heatmap_matrix))
    println(sizeof(heatmap_matrix))

    # Create a NamedArray with row and column names using sorted_keys
    named_heatmap = NamedArray(heatmap_matrix, (sorted_keys, sorted_keys), ("Rows", "Columns"))
    println(typeof(named_heatmap))
    println(sizeof(named_heatmap))
    return named_heatmap  # Return the named heatmap matrix
end


# Function to generate a heatmap with hierarchical clustering from a NamedArray
function plot_heatmap_with_clustering(named_matrix::NamedArray{Float64})
    
    # Step 1: Extract the matrix values and names from the NamedArray
    matrix_values = Matrix(named_matrix)  # Get the underlying matrix values
    row_names = names(named_matrix, 1)  # Get the row names
    col_names = names(named_matrix, 2)  # Get the column names

    # Step 2: Perform hierarchical clustering on rows and columns
    row_dist = pairwise(Euclidean(), matrix_values, dims = 1)  # Row distance matrix
    
    row_clust = hclust(row_dist, linkage=:single)  # Single linkage clustering for rows
    
    # Step 3: Reorder the matrix based on clustering
    row_order = sortperm(row_clust.order)
    col_order = sortperm(row_clust.order)
    
    ordered_matrix = matrix_values[row_order, ]
    ordered_row_names = row_names[row_order]

    # Step 4: Plot the heatmap with reordered matrix and row/column labels
    heatmap(
        ordered_matrix,
        xticks=(1:length(col_names), col_names),
        yticks=(1:length(ordered_row_names), ordered_row_names),
        xlabel="Columns",
        ylabel="Rows",
        title="Heatmap with Hierarchical Clustering",
        color=:inferno,
        aspect_ratio=:equal,
        c=:inferno,  # Inferno color map
        colorbar=true,
        yflip=true,  # Flip the y-axis to align with typical heatmap display
        size=(800, 600)
    )

    # Optionally, plot the dendrogram for rows and columns
    plot(
        plot(row_clust, xlabel="Row Dendrogram", ylabel="Distance", label="", title="Row Clustering"),
        layout=grid(2,1, heights=[0.3, 0.7])
    )
end




function plot_heatmap(heatmap_matrix::NamedArray)
    row_names = names(heatmap_matrix, 1)
    col_names = names(heatmap_matrix, 2)
    plot = heatmap(heatmap_matrix, xticks=(1:1:length(col_names), col_names), 
        yticks=(1:1:length(row_names), row_names), yflip = true,
        xlabel="Key 1", ylabel="Key 2", title="Most Significant Z-Scores Heatmap", 
        color=:inferno, xrotation=45, colorbar=true, colorbar_width=15
        )

    display(plot)
end


function zscore_heatmap(summaries::Dict{Tuple{String, String}, NamedTuple})
    most_sig_zs = extract_most_sig_z(summaries)
    heatmap_matrix = create_matrix(most_sig_zs)
    #plot_heatmap(heatmap_matrix, unique_keys)
    return heatmap_matrix
end


function export_matrix_to_csv(summaries::Dict{Tuple{String, String}, NamedTuple}; filename::Union{String, Nothing} = nothing)
    # Get the heatmap matrix and sorted keys
    sig_z_dict = extract_most_sig_z(summaries)
    heatmap_matrix = create_matrix(sig_z_dict)

    println(names(heatmap_matrix))
    # Convert the matrix to a DataFrame
    df = DataFrame(heatmap_matrix, :auto)  # Automatically generate column names (for now)

    # Rename columns with sorted keys (ensure length matches the number of columns)
    # The number of columns in `heatmap_matrix` must match the number of `sorted_keys`
    new_col_names = sorted_keys  # Concatenate 'Key1' for the row labels with sorted keys
 
    # Check if the length of new_col_names matches the number of columns in df
    if length(new_col_names) != ncol(df)
        error("The number of column names does not match the number of columns in the DataFrame.")
    end

    rename!(df, new_col_names)  # Use names! to rename all columns

    # Insert the sorted keys as the first column for row labels
    insertcols!(df, 1, :Names => sorted_keys)

    if isnothing(filename)
        filename = "zscores.csv"
    end
    # Save the DataFrame as a CSV file
    CSV.write(filename, df)
end


function export_matrix_to_csv(summaries::Dict{Tuple{String, String}, NamedTuple}; filename::Union{String, Nothing} = nothing)
    # Extract the significant Z-scores dictionary and create the heatmap matrix
    sig_z_dict = extract_most_sig_z(summaries)
    heatmap_matrix = create_matrix(sig_z_dict)  # This returns a NamedArray

    # Convert NamedArray to a regular Matrix and get row/column names
    matrix_data = Matrix(heatmap_matrix)  # Convert NamedArray to Matrix
    row_names = names(heatmap_matrix, 1)
    col_names = names(heatmap_matrix, 2)

    # Create DataFrame from the matrix and set column names
    df = DataFrame(matrix_data, Symbol.(col_names))  # Assign column names from col_names

    # Insert row names as the first column
    insertcols!(df, 1, :RowNames => row_names)

    # Define a filename if not provided
    if isnothing(filename)
        filename = "zscores.csv"
    end

    # Save the DataFrame to a CSV file
    CSV.write(filename, df)
end


function transform_to_dist_matrix(matrix::NamedArray)

    row_names = names(matrix, 1)
    col_names = names(matrix, 2)
    # Find the minimum and maximum finite values
    finite_vals = filter(!isinf, matrix[:])
    min_finite = minimum(finite_vals)
    max_finite = maximum(finite_vals)

    # Replace -Inf with min_finite and Inf with max_finite
    matrix[matrix .== -Inf] .= min_finite
    matrix[matrix .== Inf] .= max_finite

    # Center and scale
    matrix = 1 .- (matrix .- min_finite)./(max_finite .- min_finite)
    matrix = round.(matrix, digits=5)
    named_transformed_matrix = NamedArray(matrix, (row_names, col_names), ("Rows", "Columns"))

    return named_transformed_matrix
end


function transform_and_export_matrix_to_csv(matrix::Matrix{Float64}, row_names::Vector{String}, col_names::Vector{String}; filename::Union{String, Nothing} = nothing)
    # Step 1: Transform the matrix to replace Inf and NaN
    new_matrix = replace_infs_and_nans_with_extremes!(copy(matrix))  # Use a copy to avoid modifying the original

    # Step 2: Perform hierarchical clustering on the transformed matrix
    row_dist = pairwise(Euclidean(), new_matrix, dims = 1)
    row_cluster = hclust(row_dist, linkage = :complete)

    # Step 3: Reorder the rows of the original matrix and row names based on clustering order
    row_order = sortperm(row_cluster.order)
    ordered_matrix = matrix[row_order, :]  # Reorder rows in the original matrix
    ordered_row_names = row_names[row_order]  # Reorder row names based on row clustering order

    # Step 4: Convert to DataFrame with column names, and add reordered row names
    df = DataFrame(ordered_matrix, Symbol.(col_names))  # Set column names
    insertcols!(df, 1, :RowNames => ordered_row_names)  # Add row names as the first column

    # Step 5: Define a filename if not provided
    if isnothing(filename)
        filename = "transformed_ordered_matrix.csv"
    end

    # Step 6: Export to CSV
    CSV.write(filename, df)
end



function create_effect_size_csv(
    dict::Dict{Tuple{String, String}, Vector{Float64}}, 
    T::Vector{Int64}, 
    index::Int, 
    output_file::String
)
    # Extract unique row and column names
    rows = sort(unique(key[1] for key in keys(dict)))
    cols = sort(unique(key[2] for key in keys(dict)))
    
    # Initialize an empty 2D array for floats
    array = fill(0.0, length(rows), length(cols))
    
    # Populate the array with values from the dictionary
    for ((row, col), value) in dict
        row_idx = findfirst(x -> x == row, rows)
        col_idx = findfirst(x -> x == col, cols)
        array[row_idx, col_idx] = value[index]  # Access the specified index
    end
    
    # Subtract expected value from each element using the indexed T
    array .= array .- 2 * T[index]

    # Replace negative values with 0
    array .= ifelse.(array .>= 0, array, 0.0)

    # Convert to DataFrame for CSV export
    df = DataFrame()
    df[!, "Row"] = rows  # Add row names as a column
    for (j, col) in enumerate(cols)
        df[!, col] = array[:, j]
    end
    
    # Save to CSV
    CSV.write(output_file, df)
    
    return array  # Return the final 2D array
end



end # end of module