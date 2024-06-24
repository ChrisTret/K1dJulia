module DataWrangling

export process_genome_data, process_filtered_genome_data, process_genome_data_ungrouped

using CSV
using DataFrames
using Statistics


# Personal notes
#
# (string, string, string, string) two column median from path with groups
# (string, string, string) one column from path with key name
# (dataframe, string, string, string) two column median from path with groups
# (string, string, string, string, float) two column median from path with groups and filtering


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
        geno_dict[family] = sort(geno_centers)
    end
    
    return geno_dict
end

"""
    process_gennome_data(file_path::String, column_name::String, dict_key::String) -> Dict{String, Vector}

Reads a TSV file and extracts a single column of data, storing the values in a dictionary under a specified key.

# Arguments
- `file_path::String`: The path to the TSV file containing the genomic or other data.
- `column_name::String`: The name of the column from which to extract data.
- `dict_key::String`: The name to use as the key in the returned dictionary.

# Returns
- `Dict{String, Vector}`: A dictionary containing a single key-value pair, where the key is `dict_key` and the value is a vector of data extracted from `column_name`.
"""
function process_genome_data(file_path::String, column_name::String, dict_key::String)::Dict{String, Vector}
    # Read the TSV file into a DataFrame
    data = CSV.read(file_path, DataFrame; header = 1)

    # Extract the column data
    column_data = data[!, column_name]

    # Create a dictionary to store the results
    result_dict = Dict{String, Vector}()
    result_dict[dict_key] = sort(column_data)

    return result_dict
end


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
function process_genome_data(data::DataFrame, cluster_group::String, start_colname::String, end_colname::String)::Dict{String, Vector{Int64}}
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
        geno_dict[family] = sort(geno_centers)
    end
    
    return geno_dict
end



"""
    process_genome_data(file_path::String, cluster_group::String, start_colname::String, end_colname::String, percentile::Float64) -> Dict{String, Vector{Int64}}

Processes genomic data by filtering to keep only the longest `percentile` entries for each `cluster_group`.

# Arguments
- `file_path::String`: Path to the CSV file containing the genomic data.
- `cluster_group::String`: Column name to group by (e.g., 'repFamily').
- `start_colname::String`: Column name indicating the start of an element.
- `end_colname::String`: Column name indicating the end of an element.
- `percentile::Float64`: Percentile value to filter by, within each group.

# Returns
- `Dict{String, Vector{Int64}}`: Dictionary with each group as keys and vectors of genomic centers of the longest entries as values.
"""
function process_genome_data(file_path::String, cluster_group::String, start_colname::String, end_colname::String, percentile::Float64)::Dict{String, Vector{Int64}}
    data = CSV.read(file_path, DataFrame; header = 1)
    # Calculate length for each entry
    data[!, "Length"] = data[:, end_colname] - data[:, start_colname]

    # Initialize dictionary to store results
    filtered_dict = Dict{String, Vector{Int64}}()

    # Group by `cluster_group` and filter within each group
    grouped_data = groupby(data, cluster_group)
    for group in grouped_data
        # Calculate cutoff within each group
        cutoff = quantile(group[!, "Length"], 1 - percentile)
        filtered_group = filter(row -> row["Length"] >= cutoff, group)
        
        # Calculate centers for filtered data
        filtered_group[!, "Center"] = round.(Int, (filtered_group[:, start_colname] .+ filtered_group[:, end_colname]) ./ 2)
        
        # Store the centers in the dictionary
        filtered_dict[group[1, cluster_group]] = sort(filtered_group[!, "Center"])
    end
    
    return filtered_dict
end

"""
    process_genome_centers(file_path::String, start_colname::String, end_colname::String, dict_key::String) -> Dict{String, Vector{Int64}}

Reads a TSV file containing genomic data, calculates the `genoCenter` for each entry, and stores these `genoCenter` values in a dictionary under a specified key.

# Arguments
- `file_path::String`: The path to the TSV file containing the genomic data.
- `start_colname::String`: The name of the column indicating the start of an element.
- `end_colname::String`: The name of the column indicating the end of an element.
- `dict_key::String`: The key under which the `genoCenter` values will be stored in the returned dictionary.

# Returns
- `Dict{String, Vector{Int64}}`: A dictionary containing a single key-value pair, where the key is provided by `dict_key` and the value is a vector of `genoCenter` values calculated from the data.
"""
function process_genome_data_ungrouped(file_path::String, start_colname::String, end_colname::String, dict_key::String)::Dict{String, Vector{Int64}}
    # Read the TSV file into a DataFrame
    data = CSV.read(file_path, DataFrame; header = 1)

    # Calculate genoCenter and convert to integers
    centers = sort(round.(Int, (data[:, start_colname] .+ data[:, end_colname]) ./ 2))

    # Create a dictionary to store the results
    result_dict = Dict{String, Vector{Int64}}(dict_key => centers)

    return result_dict
end


end # end of module