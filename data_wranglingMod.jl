module DataWrangling

export process_genome_data, process_filtered_genome_data

using CSV
using DataFrames
using Statistics

"""
    process_genome_data(file_path::String) -> Dict{String, Vector{Float64}}

Reads a TSV file containing genomic data, calculates the `genoCenter` for each entry, and organizes the `genoCenter` values by `repFamily`.

# Arguments
- `file_path::String`: The path to the TSV file containing the genomic data.
- `cluster_group::String`: The elements to be grouped over (eg. Alu, L1, etc.)
- `start_colname::String`: The name of the column indicating the start of an element
- `end_colname::String`: The name of the column indicating the end of an element

# Returns
- `Dict{String, Vector{Int64}}`: A dictionary where each key is a unique `repFamily` and each value is a vector of `genoCenter` values associated with that `repFamily`.
"""
function process_genome_data(file_path::String, cluster_group::String, start_colname::String, end_colname::String)::Dict{String, Vector{Int64}}
    # Read the TSV file into a DataFrame
    data = CSV.read(file_path, DataFrame; header = 1)

    # Convert genoCenter to integers
    data[!, "Center"] = round.(Int, (data[:, start_colname] .+ data[:, end_colname]) ./ 2)
    # Initialize an empty dictionary to store the results
    geno_dict = Dict{String, Vector{Int64}}()
    
    # Group the data by repFamily
    grouped_data = groupby(data, cluster_group)
    
    # Process each group
    for g in grouped_data
        # Get the repFamily name (all entries in g have the same repFamily)
        family = g[1, cluster_group]
        # Get the genoCenter values for this group
        geno_centers = g[:, "Center"]
        # Store the genoCenters in the dictionary
        geno_dict[family] = geno_centers
    end
    
    return geno_dict
end


"""
    process_filtered_genome_data(file_path::String, cluster_group::String, start_colname::String, end_colname::String, percentile::Float64) -> Dict{String, Vector{Int64}}

Read a CSV file containing genomic data, filter entries based on the length percentile specified by the user, and compute the genomic center for each entry, organized by the specified `cluster_group`.

# Arguments
- `file_path::String`: The path to the CSV file containing the genomic data.
- `cluster_group::String`: The name of the column by which the data should be grouped (e.g., 'repFamily').
- `start_colname::String`: The name of the column indicating the start of an element.
- `end_colname::String`: The name of the column indicating the end of an element.
- `percentile::Float64`: The percentile value used to filter data based on the length of elements. The function filters out the shortest elements, keeping only those in the top `1 - percentile` fraction.

# Returns
- `Dict{String, Vector{Int64}}`: A dictionary where each key corresponds to a unique value from the `cluster_group` column, and each value is a vector of integers representing the genomic center of the filtered entries for that group.
"""
function process_filtered_genome_data(file_path::String, cluster_group::String, start_colname::String, end_colname::String, percentile::Float64)::Dict{String, Vector{Int64}}
    data = CSV.read(file_path, DataFrame; header = 1)
    # Calculate lengths and filter based on the specified percentile
    data[!, "Length"] = data[:, end_colname] - data[:, start_colname]
    cutoff = quantile(data[!, "Length"], 1 - percentile)
    filtered_data = filter(row -> row["Length"] >= cutoff, data)
    # Calculate centers of the filtered data
    filtered_data[!, "Center"] = round.(Int, (filtered_data[:, start_colname] .+ filtered_data[:, end_colname]) ./ 2)
    filtered_dict = Dict{String, Vector{Int64}}()
    grouped_data = groupby(filtered_data, cluster_group)
    for g in grouped_data
        family = g[1, cluster_group]
        geno_centers = g[:, "Center"]
        filtered_dict[family] = geno_centers
    end
    return filtered_dict
end


end # end of module