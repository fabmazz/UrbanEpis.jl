module UrbanEpis


using Random

using DataFrames
using StatsBase

using Graphs
using StaticArrays: SVector

export draw_people_same_tile, extract_tile_data, SMIRModel, SMIRStatesCounter, filter_tidx_by_pop, build_cells_neighbors_graph

include("city.jl")
include("plot_utils.jl")
include("dfutils.jl")
include("epidemic.jl")
include("epipars.jl")
include("stats.jl")
include("stats_seir.jl")

export df_fields_op, rate_to_prob
export find_fraction_infected, draw_misbeh_hom,draw_misbeh_hom, run_epidemics_city, run_epidemics_city_parallel

export find_avg_attack_rates_allsims, tiles_misinfect_stats_full, calc_curves_misinf

export calc_frac_infected_tile_SMIR, save_trace_inf_misinf, count_SMIR_state

export calc_tiles_timeinf_infector, calc_tiles_inf_history, calc_pmisinf_from_pvhe, calc_flow_infection_tiles, join_all_flow_infect

export count_states_by_tileid, dictForArr, count_people_active_tile, bin_counts_tileids

end
