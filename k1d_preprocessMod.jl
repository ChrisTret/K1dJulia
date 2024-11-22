module K1dPreprocess

using CSV, Printf, Plots, DataFrames, Measures

export display_counts, filter_data_by_occurrences, count_non_empty_regions, print_regions_with_most_entries, occurrences_plot, calculate_element_lengths
"""
    display_counts(data_dict::Dict{String, Dict{String, Vector{Int}}})

Calculates and displays the total counts of values for each name in the provided dictionary of dictionaries. The output is sorted in descending order of total counts.

# Arguments
- `data_dict::Dict{String, Dict{String, Vector{Int}}}`: A dictionary where each key is a name (e.g., gene name), and each value is another dictionary. The inner dictionary maps chromosome names to vectors of integers.

# Returns
- `Nothing`: This function prints the results directly to the console and does not return any value.
"""
function display_counts(data_dict::Dict{String, Dict{String, Vector{Int}}})
    total_counts = Dict{String, Int}()
    
    for (name, chrom_dict) in data_dict
        total_counts[name] = sum(length(values) for values in values(chrom_dict))
    end
    
    sorted_names = sort(collect(total_counts), by = x -> x[2], rev=true)
    
    # Display the results
    println(@sprintf("%-20s %-10s", "Name", "Total"))
    println("=" ^ 30)
    for (name, total) in sorted_names
        println(@sprintf("%-20s %-10d", name, total))
    end
end


"""
    filter_data_by_occurrences!(data_dict::Dict{String, Dict{String, Vector{Int}}}, min_occurrences::Int) 
    -> Dict{String, Dict{String, Vector{Int}}}

Filters the provided dictionary of data in place, keeping only entries that have a total number of occurrences across all chromosomes greater than or equal to a specified minimum.

# Arguments
- `data_dict::Dict{String, Dict{String, Vector{Int}}}`: A dictionary where each key is a name (e.g., gene name), and each value is another dictionary. The inner dictionary maps chromosome names to vectors of integers.
- `min_occurrences::Int`: The minimum total number of occurrences required for a name to be kept in the dictionary.

# Returns
- `Dict{String, Dict{String, Vector{Int}}}`: The filtered dictionary containing only entries with total occurrences greater than or equal to `min_occurrences`.
"""
function filter_data_by_occurrences(
    data_dict::Dict{String, Dict{String, Vector{Int}}},
    min_occurrences::Union{Int, Nothing} = nothing,
    max_occurrences::Union{Int, Nothing} = nothing
)
    # Create a new dictionary to store the filtered data
    filtered_dict = Dict{String, Dict{String, Vector{Int}}}()
    
    for (name, chrom_dict) in data_dict
        # Calculate the total number of occurrences across all chromosomes
        total_count = sum(length(values) for values in values(chrom_dict))
        
        # Apply filtering conditions
        if (!isnothing(min_occurrences) && total_count < min_occurrences) || 
           (!isnothing(max_occurrences) && total_count > max_occurrences)
            continue
        end
        
        # If it meets the condition, add it to the filtered dictionary
        filtered_dict[name] = chrom_dict
    end
    
    # Return the filtered dictionary
    return filtered_dict
end





function count_non_empty_regions(data_dict::Dict{String, Dict{Int, Vector{Int64}}})
    # Dictionary to store the count of non-empty regions for each name
    non_empty_count = Dict{String, Int}()
    # Dictionary to store the region sizes for each name in descending order
    regions_by_size = Dict{String, Vector{Tuple{Int, Int}}}()

    for (name, regions) in data_dict
        # Find non-empty regions and their sizes
        non_empty_regions = filter(x -> !isempty(x[2]), regions)
        non_empty_count[name] = length(non_empty_regions)

        # Sort regions by size in descending order
        sorted_regions = sort(collect(non_empty_regions), by = x -> -length(x[2]))
        regions_by_size[name] = [(region[1], length(region[2])) for region in sorted_regions]
    end

    # Sort names by the number of non-empty regions in descending order
    sorted_names = sort(collect(keys(non_empty_count)), by = x -> -non_empty_count[x])

    # Print the names in the order of non-empty regions count
    for name in sorted_names
        println("$name has ", non_empty_count[name], " non-empty regions.")
    end

    return regions_by_size
end

# Example usage:
# result_dict = process_genome_data_by_region("path/to/file.tsv", name_col="Name", num_regions=5, start_col="Start", end_col="End")
# regions_data = count_non_empty_regions(result_dict)

# To print regions with the most entries for a given name
function print_regions_with_most_entries(regions_data::Dict{String, Vector{Tuple{Int, Int}}}, name::String)
    if haskey(regions_data, name)
        println("Regions for $name with most entries:")
        for (region, count) in regions_data[name]
            println("Region $region: $count entries")
        end
    else
        println("No data found for name: $name")
    end
end


function occurrences_plot(genome_data_dict::Dict{String, Dict{String, Vector{Int64}}}, title ::Union{String, Nothing} = nothing)
    # Prepare storage for the key names and their respective counts
    keys_list = String[]
    counts_list = Int[]

    # Loop over the outer dictionary (names or chroms)
    for (outer_key, inner_dict) in genome_data_dict
        total_count = 0
        
        # Loop over the inner dictionary (chromosomes) to count the total occurrences
        for (chrom_key, values) in inner_dict
            total_count += length(values)  # Add the number of occurrences for this chrom
        end
        
        # Store the outer key and the corresponding total count
        push!(keys_list, outer_key)
        push!(counts_list, total_count)
    end

    # Combine the keys and counts into a tuple and sort them by counts
    sorted_data = sort(collect(zip(keys_list, counts_list)), by = x -> x[2])

    # Unzip the sorted data
    sorted_keys = [x[1] for x in sorted_data]
    sorted_counts = [x[2] for x in sorted_data]

    if isnothing(title)
        title = "Occurrences per Key"
    end
    # Create a bar plot with logarithmic y-axis
    bar_plot = bar(sorted_keys, sorted_counts, yscale=:log10, xlabel="Classes", ylabel="Occurrences (log scale)", 
        title=title, legend=false, xrotation=45, margin = 5mm, fmt = :png)

    display(bar_plot)
end


function calculate_element_lengths(
    file_path::String; 
    name_col::Union{String, Nothing}=nothing, 
    chrom_col::String, 
    start_col::String, 
    end_col::String
)::Dict{String, Vector{Int64}}
    
    # Read the TSV file into a DataFrame
    data = CSV.read(file_path, DataFrame; header = 1)
    result_dict = Dict{String, Vector{Int64}}()

    # Calculate the length of each element as end - start
    data[!, "Length"] = data[:, end_col] .- data[:, start_col]

    if name_col !== nothing
        # Group by name and chromosome
        grouped_data = groupby(data, [name_col, chrom_col])
        for g in grouped_data
            name = g[1, name_col]
            lengths = g[:, "Length"]

            if !haskey(result_dict, name)
                result_dict[name] = Vector{Int64}()
            end
            append!(result_dict[name], lengths)
        end
    else
        # Group by chromosome only
        grouped_data = groupby(data, chrom_col)
        for g in grouped_data
            chrom = g[1, chrom_col]
            lengths = g[:, "Length"]

            if !haskey(result_dict, chrom)
                result_dict[chrom] = Vector{Int64}()
            end
            append!(result_dict[chrom], lengths)
        end
    end

    return result_dict
end


# function occurrences_plot(genome_data_dict::Dict{String, Dict{String, Vector{Int64}}})
#     # Prepare storage for the key names and their respective counts
#     keys_list = String[]
#     counts_list = Int[]

#     # Loop over the outer dictionary (names or chroms)
#     for (outer_key, inner_dict) in genome_data_dict
#         total_count = 0
        
#         # Loop over the inner dictionary (chromosomes) to count the total occurrences
#         for (chrom_key, values) in inner_dict
#             total_count += length(values)  # Add the number of occurrences for this chrom
#         end
        
#         # Store the outer key and the corresponding total count
#         push!(keys_list, outer_key)
#         push!(counts_list, total_count)
#     end

#     # Combine the keys and counts into a tuple and sort them by counts
#     sorted_data = sort(collect(zip(keys_list, counts_list)), by = x -> x[2])

#     # Unzip the sorted data
#     sorted_keys = [x[1] for x in sorted_data]
#     sorted_counts = [x[2] for x in sorted_data]

#     # Create a bar plot with a logarithmic y-axis using PlotlyJS
#     bar_plot = Plot(
#         data=[bar(x=sorted_keys, y=sorted_counts, text=sorted_counts, textposition="auto")],
#         Layout(
#             title="Occurrences per Key",
#             xaxis_title="Keys",
#             yaxis_title="Occurrences (log scale)",
#             yaxis=attr(type="log"),
#             bargap=0.2,
#             xaxis=attr(tickangle=-45)
#         )
#     )

#     display(bar_plot)
# end


function length_for_key(
    data::Dict{String, Dict{String, Vector{Int64}}}, key::String
) ::Int
    # Check if the key exists in the dictionary
    if haskey(data, key)
        # Calculate the sum of lengths of vectors in the inner dictionary for the specified key
        return sum(length(vec) for vec in values(data[key]))
    else
        error("Key not found in dictionary")
    end
end

end # end of module