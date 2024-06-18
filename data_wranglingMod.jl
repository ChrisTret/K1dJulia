module DataWrangling

export process_genome_data

using CSV
using DataFrames

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

end # end of module