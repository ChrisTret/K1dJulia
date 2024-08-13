module K1dPreprocess

using CSV, Printf

export display_counts, filter_data_by_occurrences!

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