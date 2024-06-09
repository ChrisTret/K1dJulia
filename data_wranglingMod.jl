module DataWrangling

export process_genome_data

using CSV
using DataFrames

"""
    process_genome_data(file_path::String) -> Dict{String, Vector{Float64}}

Reads a TSV file containing genomic data, calculates the `genoCenter` for each entry, and organizes the `genoCenter` values by `repFamily`.

# Arguments
- `file_path::String`: The path to the TSV file containing the genomic data.

# Returns
- `Dict{String, Vector{Int64}}`: A dictionary where each key is a unique `repFamily` and each value is a vector of `genoCenter` values associated with that `repFamily`.
"""
function process_genome_data(file_path::String)::Dict{String, Vector{Int64}}
    # Read the TSV file into a DataFrame
    data = CSV.read(file_path, DataFrame; delim='\t')

     # Convert genoCenter to integers
    data[!, :genoCenter] = round.(Int, (data[:, :genoStart] .+ data[:, :genoEnd]) ./ 2)
    # Initialize an empty dictionary to store the results
    geno_dict = Dict{String, Vector{Int64}}()
    
    # Group the data by repFamily
    grouped_data = groupby(data, :repFamily)
    
    # Process each group
    for g in grouped_data
        # Get the repFamily name (all entries in g have the same repFamily)
        family = g[1, :repFamily]
        # Get the genoCenter values for this group
        geno_centers = g[:, :genoCenter]
        # Store the genoCenters in the dictionary
        geno_dict[family] = geno_centers
    end
    
    return geno_dict
end

end # end of module