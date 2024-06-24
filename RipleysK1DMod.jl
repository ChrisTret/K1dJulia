module RipleysK1D

using CSV, DataFrames, Statistics, Plots

include("data_wranglingMod.jl")
include("k1d_functionsMod.jl")
include("k1d_speed_testsMod.jl")
include("k_plotsMod.jl")
include("k1d_conf_intMod.jl")
include("significance_measuresMod.jl")

using .DataWrangling, .K1dFun, .K1dSpeedTests, .K1dPlots, .K1dConfInt, .SigMeasures

export process_genome_data, merge_dictionaries
export k1d_univ, k1d_biv, k1d_all_comparisons 
export speed_comparison_univ 
export k_plot, plot_top_n, plot_by_key, plot_pair
export bootstrap_variance_k1d 
export max_percent_increase

end # end of module

