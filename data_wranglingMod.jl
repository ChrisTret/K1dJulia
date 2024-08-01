module DataWrangling

export process_genome_data, merge_dictionaries, process_genome_data_by_chromosome, process_genome_data_by_region, process_genome_data_by_chromosome_by_region

using CSV
using DataFrames
using Statistics

"""
    process_genome_data(file_path::String; cluster_group::Union{String, Nothing}=nothing, start_col::Union{String, Nothing}=nothing, end_col::Union{String, Nothing}=nothing, dict_key::Union{String, Nothing}=nothing, percentile::Union{Float64, Nothing}=nothing) -> Dict{String, Vector{Int64}}

Processes genomic data from a TSV file according to various scenarios.

# Arguments
- `file_path::String`: The path to the TSV file containing the genomic data.
- `cluster_group::Union{String, Nothing}`: The column name to group by (e.g., 'repFamily').
- `start_col::Union{String, Nothing}`: The name of the column indicating the start of an element.
- `end_col::Union{String, Nothing}`: The name of the column indicating the end of an element.
- `dict_key::Union{String, Nothing}`: The name to use as the key in the returned dictionary.

# Returns
- `Dict{String, Vector{Int64}}`: A dictionary where each key is a group name or `dict_key`, and each value is a vector of calculated centers or column values.
"""
function process_genome_data(file_path::String; cluster_group::Union{String, Nothing}=nothing, start_col::Union{String, Nothing}=nothing, end_col::Union{String, Nothing}=nothing, dict_key::Union{String, Nothing}=nothing)::Dict{String, Vector{Int64}}
    # Read the TSV file into a DataFrame
    data = CSV.read(file_path, DataFrame; header = 1)
    result_dict = Dict{String, Vector{Int64}}()

    if start_col !== nothing && end_col !== nothing
        # Calculate genoCenter
        data[!, "Center"] = round.(Int, (data[:, start_col] .+ data[:, end_col]) ./ 2)

        if cluster_group !== nothing
            # Group the data by cluster_group
            grouped_data = groupby(data, cluster_group)
            for g in grouped_data
                family = g[1, cluster_group]
                geno_centers = g[:, "Center"]
                result_dict[family] = sort(geno_centers)
            end
        else
            # No grouping, use dict_key or default key
            key = dict_key === nothing ? "all_data" : dict_key
            result_dict[key] = sort(data[!, "Center"])
        end
    elseif start_col !== nothing
        # Process data with only start column

        if cluster_group !== nothing
            # Group the data by cluster_group
            grouped_data = groupby(data, cluster_group)
            for g in grouped_data
                family = g[1, cluster_group]
                geno_starts = g[:, start_col]
                result_dict[family] = sort(geno_starts)
            end
        else
            # No grouping, use dict_key or default key
            key = dict_key === nothing ? "all_data" : dict_key
            result_dict[key] = sort(data[!, start_col])
        end
    else
        # Reject attempt
        println("Invalid function inputs for process_genome_data")
        result_dict = nothing
    end

    return result_dict
end


"""
    process_genome_data_by_chromosome(file_path::String, name_col::String, chrom_col::String, start_col::Union{String, Nothing}=nothing, end_col::Union{String, Nothing}=nothing) -> Dict{String, Dict{String, Vector{Int64}}}

Processes genomic data from a TSV file and organizes it by name and chromosome.

# Arguments
- `file_path::String`: The path to the TSV file containing the genomic data.
- `name_col::String`: The column name for the genomic element (e.g., 'repName').
- `chrom_col::String`: The column name for the chromosome (e.g., 'chrom').
- `start_col::Union{String, Nothing}`: The name of the column indicating the start of an element.
- `end_col::Union{String, Nothing}`: The name of the column indicating the end of an element.

# Returns
- `Dict{String, Dict{String, Vector{Int64}}}`: A nested dictionary where the first key is the name, and the second key is the chromosome, with each value being a vector of calculated centers or column values.
"""
function process_genome_data_by_chromosome(file_path::String; name_col::Union{String, Nothing}=nothing, chrom_col::String, start_col::Union{String, Nothing}=nothing, end_col::Union{String, Nothing}=nothing, dict_key::Union{String, Nothing}=nothing)::Dict{String, Dict{String, Vector{Int64}}}
    # Read the TSV file into a DataFrame
    data = CSV.read(file_path, DataFrame; header = 1)
    result_dict = Dict{String, Dict{String, Vector{Int64}}}()

    if start_col !== nothing && end_col !== nothing
        # Calculate genoCenter
        data[!, "Center"] = round.(Int, (data[:, start_col] .+ data[:, end_col]) ./ 2)
        data_col = "Center"
    elseif start_col !== nothing
        # Use start column
        data_col = start_col
    else
        # Reject attempt
        println("Invalid function inputs for process_genome_data_by_chromosome")
        return nothing
    end

    if name_col !== nothing
        # Group by name and chromosome
        grouped_data = groupby(data, [name_col, chrom_col])
        for g in grouped_data
            name = g[1, name_col]
            chrom = g[1, chrom_col]
            values = sort(g[:, data_col])

            if !haskey(result_dict, name)
                result_dict[name] = Dict{String, Vector{Int64}}()
            end
            result_dict[name][chrom] = values
        end
    else
        # Group by chromosome only
        grouped_data = groupby(data, chrom_col)
        key = dict_key === nothing ? "all_data" : dict_key

        for g in grouped_data
            chrom = g[1, chrom_col]
            values = sort(g[:, data_col])

            if !haskey(result_dict, key)
                result_dict[key] = Dict{String, Vector{Int64}}()
            end
            result_dict[key][chrom] = values
        end
    end

    return result_dict
end



function process_genome_data_by_region(file_path::String; name_col::String, num_regions::Int, chrom_col::Union{String, Nothing}=nothing, start_col::Union{String, Nothing}=nothing, end_col::Union{String, Nothing}=nothing)::Dict{String, Dict{Int, Vector{Int64}}}
    # Read the TSV file into a DataFrame
    data = CSV.read(file_path, DataFrame; header = 1)
    result_dict = Dict{String, Dict{Int, Vector{Int64}}}()

    if start_col !== nothing && end_col !== nothing
        # Calculate genoCenter
        data[!, "Center"] = round.(Int, (data[:, start_col] .+ data[:, end_col]) ./ 2)
        data_col = "Center"
    elseif start_col !== nothing
        # Use start column
        data_col = start_col
    else
        # Reject attempt
        println("Invalid function inputs for process_genome_data_by_region")
        return nothing
    end

    if chrom_col !== nothing
        # Group by name and chromosome
        grouped_data = groupby(data, [name_col, chrom_col])
        for g in grouped_data
            name = g[1, name_col]
            chrom = g[1, chrom_col]
            values = sort(g[:, data_col])

            # Determine the min and max of the data column
            min_val = minimum(values)
            max_val = maximum(values)
            region_size = ceil(Int, (max_val - min_val + 1) / num_regions)

            # Initialize the nested dictionaries if they don't exist
            if !haskey(result_dict, "$name:$chrom")
                result_dict["$name:$chrom"] = Dict{Int, Vector{Int64}}()
            end

            # Allocate regions
            for i in 1:num_regions
                result_dict["$name:$chrom"][i] = Int64[]
            end

            # Assign values to regions
            for val in values
                region = min(num_regions, ceil(Int, (val - min_val + 1) / region_size))
                push!(result_dict["$name:$chrom"][region], val)
            end
        end
    else
        # Group by name only
        grouped_data = groupby(data, name_col)
        for g in grouped_data
            name = g[1, name_col]
            values = sort(g[:, data_col])

            # Determine the min and max of the data column
            min_val = minimum(values)
            max_val = maximum(values)
            region_size = ceil(Int, (max_val - min_val + 1) / num_regions)

            # Initialize the nested dictionary if it doesn't exist
            if !haskey(result_dict, name)
                result_dict[name] = Dict{Int, Vector{Int64}}()
            end

            # Allocate regions
            for i in 1:num_regions
                result_dict[name][i] = Int64[]
            end

            # Assign values to regions
            for val in values
                region = min(num_regions, ceil(Int, (val - min_val + 1) / region_size))
                push!(result_dict[name][region], val)
            end
        end
    end

    return result_dict
end

function process_genome_data_by_chromosome_by_region(file_path::String; name_col::String, num_regions::Int, chrom_col::String, start_col::Union{String, Nothing}=nothing, end_col::Union{String, Nothing}=nothing)::Dict{String, Dict{String, Dict{Int, Vector{Int64}}}} 
    # Read the TSV file into a DataFrame
    data = CSV.read(file_path, DataFrame; header = 1)
    result_dict = Dict{String, Dict{String, Dict{Int, Vector{Int64}}}}()

    if start_col !== nothing && end_col !== nothing
        # Calculate genoCenter
        data[!, "Center"] = round.(Int, (data[:, start_col] .+ data[:, end_col]) ./ 2)
        data_col = "Center"
    elseif start_col !== nothing
        # Use start column
        data_col = start_col
    else
        # Reject attempt
        println("Invalid function inputs for process_genome_data_by_chromosome_by_region")
        return nothing
    end

    # Group by name and chromosome
    grouped_data = groupby(data, [name_col, chrom_col])
    for g in grouped_data
        name = g[1, name_col]
        chrom = g[1, chrom_col]
        values = sort(g[:, data_col])

        # Determine the min and max of the data column
        min_val = minimum(values)
        max_val = maximum(values)
        region_size = ceil(Int, (max_val - min_val + 1) / num_regions)

        # Initialize the nested dictionaries if they don't exist
        if !haskey(result_dict, name)
            result_dict[name] = Dict{String, Dict{Int, Vector{Int64}}}()
        end
        if !haskey(result_dict[name], chrom)
            result_dict[name][chrom] = Dict{Int, Vector{Int64}}()
        end

        # Allocate regions
        for i in 1:num_regions
            result_dict[name][chrom][i] = Int64[]
        end

        # Assign values to regions
        for val in values
            region = min(num_regions, ceil(Int, (val - min_val + 1) / region_size))
            push!(result_dict[name][chrom][region], val)
        end
    end

    return result_dict
end


"""
    merge_dictionaries(dicts::Vector{Dict{String, Vector{Int64}}}) -> Dict{String, Vector{Int64}}

Merges multiple dictionaries into a single dictionary, handling duplicate keys by renaming them and printing a summary of duplicates.

# Arguments
- `dicts::Vector{Dict{String, Vector{Int64}}}`: A vector of dictionaries to merge.

# Returns
- `Dict{String, Vector{Int64}}`: A dictionary with merged keys and concatenated vectors.
"""
function merge_dictionaries(dicts::Vector{Dict{String, Vector{Int64}}})::Dict{String, Vector{Int64}}
    merged_dict = Dict{String, Vector{Int64}}()
    duplicate_info = DataFrame(Name=String[], NewName=String[], Dictionary=String[])
    key_counts = Dict{String, Int}()

    for (i, dict) in enumerate(dicts)
        for (key, value) in dict
            new_key = key
            if haskey(key_counts, key)
                key_counts[key] += 1
                new_key = "$(key)_duplicate_$(key_counts[key])"
                push!(duplicate_info, (Name=key, NewName=new_key, Dictionary="dict$i"))
            else
                key_counts[key] = 0
            end

            if haskey(merged_dict, new_key)
                append!(merged_dict[new_key], value)
            else
                merged_dict[new_key] = copy(value)
            end
        end
    end

      if nrow(duplicate_info) > 0
        println("Duplicate Keys Summary:")
        println(duplicate_info)
    end

    return merged_dict
end


using DataFrames

function merge_dictionaries(dicts::Vector{Dict{String, Dict{String, Vector{Int64}}}})::Dict{String, Dict{String, Vector{Int64}}}
    merged_dict = Dict{String, Dict{String, Vector{Int64}}}()
    duplicate_info = DataFrame(Name=String[], NewName=String[], Dictionary=String[])
    key_counts = Dict{String, Int}()

    for (i, dict) in enumerate(dicts)
        for (outer_key, inner_dict) in dict
            new_outer_key = outer_key
            if haskey(key_counts, outer_key)
                key_counts[outer_key] += 1
                new_outer_key = "$(outer_key)_duplicate_$(key_counts[outer_key])"
                push!(duplicate_info, (Name=outer_key, NewName=new_outer_key, Dictionary="dict$i"))
            else
                key_counts[outer_key] = 0
            end

            if !haskey(merged_dict, new_outer_key)
                merged_dict[new_outer_key] = Dict{String, Vector{Int64}}()
            end

            for (inner_key, value) in inner_dict
                if haskey(merged_dict[new_outer_key], inner_key)
                    append!(merged_dict[new_outer_key][inner_key], value)
                else
                    merged_dict[new_outer_key][inner_key] = copy(value)
                end
            end
        end
    end

    if nrow(duplicate_info) > 0
        println("Duplicate Keys Summary:")
        println(duplicate_info)
    end

    return merged_dict
end


end # end of module
