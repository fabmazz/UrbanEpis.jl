module UrbanEpis


using Random

using DataFrames
using StatsBase

using Graphs
#using StatsBase

export draw_people_same_tile, extract_tile_data

include("city.jl")
include("plot_utils.jl")
include("dfutils.jl")
include("epidemic.jl")
include("stats.jl")

export df_fields_op
export draw_misinformed, find_fraction_infected, draw_misinformed_tile, draw_misinformed_homogen, run_epidemics_city, run_epidemics_city_parallel

export find_avg_attack_rates_allsims, tiles_misinfect_stats_full

end
