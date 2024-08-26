module K1dPreprocess

using CSV, Printf

export display_counts, filter_data_by_occurrences!
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
function filter_data_by_occurrences!(data_dict::Dict{String, Dict{String, Vector{Int}}}, min_occurrences::Int)
    # Create a new dictionary to store the filtered data
    filtered_dict = Dict{String, Dict{String, Vector{Int}}}()
    
    for (name, chrom_dict) in data_dict
        # Calculate the total number of occurrences across all chromosomes
        total_count = sum(length(values) for values in values(chrom_dict))
        
        # If the total count is greater than or equal to the minimum, keep it
        if total_count >= min_occurrences
            filtered_dict[name] = chrom_dict
        end
    end
    
    # Update the original dictionary with the filtered data
    empty!(data_dict)
    merge!(data_dict, filtered_dict)
    
    return data_dict
end

end # end of module