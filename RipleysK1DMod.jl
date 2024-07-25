module RipleysK1D

using CSV, DataFrames, Statistics, Plots

include("data_wranglingMod.jl")
include("k1d_functionsMod.jl")
include("k1d_speed_testsMod.jl")
include("k_plotsMod.jl")
include("k1d_conf_intMod.jl")
include("significance_measuresMod.jl")
include("k1d_monte_carlo.jl")

using .DataWrangling, .K1dFun, .K1dSpeedTests, .K1dPlots, .K1dConfInt, .SigMeasures, .K1dMonteCarlo

export process_genome_data, merge_dictionaries, process_genome_data_by_chromosome, process_genome_data_by_region, process_genome_data_by_chromosome_by_region
export k1d_univ, k1d_biv, k1d_all_comparisons, k1d_all_comparisons_indiv_chrom,  k1d_mean_across_chromosomes, k1d_mean_across_regions
export speed_comparison_univ 
export k_plot, l_plot, plot_top_n, plot_by_key, plot_pair
export bootstrap_variance_k1d, bootstrap_variance_k1d_all_comparisons, approx_k1d_mean_variance_chromosome
export max_percent_increase
export monte_carlo_sim, monte_carlo_test, extract_mcecb

end # end of module