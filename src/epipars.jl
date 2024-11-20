export HesiTileData, shuffle_hesit_tile!,calc_tile_pm_hes_lin

rate_to_prob(b::Real,t::Real) = 1 - exp(-b*t)

struct CitySimPars1
    R0::AbstractFloat
    ts::AbstractFloat
    a::Real
    b::Real
end

struct HesiTileData{I<:Integer, F<: AbstractFloat}
    #fixed_data::TileData{I}
    tile_id::Vector{I}
    tile_pop::Vector{I}
    tile_h::Vector{F}
end

CitySimPars = CitySimPars1

function calc_pmisinf_from_pvhe(datatile::TileData, tiles_vhes_to_filt::AbstractDataFrame, total_people::Integer, desired_mean_pmis::AbstractFloat; col_vaxhes="pvhe")
    df = filter("tile_id" => x -> in(x, datatile.tiles_idcs), tiles_vhes_to_filt)
    sort!(df, "tile_id")
    
    x = replace(df[!,col_vaxhes], missing => NaN)
    @assert any(isnan.(x)) == false
    ## compute cumulative distribution on selected tiles
    f_s = ecdf(x)
    upm = f_s(x) ### These are from 0 to 1
    ## calculate average unscaled_pm
    AV_upm = sum(datatile.tiles_pop[tid]*y for (tid,y) in zip(df[!,"tile_id"], upm) ) / total_people
    Mscale = desired_mean_pmis / AV_upm
    pmis_tile_order_vhes =  upm * Mscale
    return DataFrame(tile_id=df[!,"tile_id"], p_mis=pmis_tile_order_vhes)
end

function max_global_pmis_from_pvhe(datatile::TileData, tiles_vhes_to_filt::AbstractDataFrame, total_people::Integer; col_vaxhes="pvhe")
    df = filter("tile_id" => x -> in(x, datatile.tiles_idcs), tiles_vhes_to_filt)
    sort!(df, "tile_id")
    
    x = replace(df[!,col_vaxhes], missing => NaN)
    @assert any(isnan.(x)) == false
    ## compute cumulative distribution on selected tiles
    f_s = ecdf(x)
    upm = f_s(x) ### These are from 0 to 1
    ## calculate average unscaled_pm
    AV_upm = sum(datatile.tiles_pop[tid]*y for (tid,y) in zip(df[!,"tile_id"], upm) ) / total_people

    return AV_upm
end

function HesiTileData(tiledata::TileData,  tiles_vhes_to_filt::AbstractDataFrame;  col_vaxhes="pvhe")
    df = filter("tile_id" => x -> in(x, tiledata.tiles_idcs), tiles_vhes_to_filt)
    sort!(df, "tile_id")

    x = replace(df[!,col_vaxhes], missing => NaN)
    @assert any(isnan.(x)) == false
    people = collect(tiledata.tiles_pop[i] for i in df[!,:tile_id])
    idcs = convert(Int, df[!,:tile_id])
    @assert all(idcs .== tiledata.tiles_idcs)
    upm = (x.-minimum(x))./(maximum(x)-minimum(x))
    HesiTileData(idcs,people, upm)
end

function shuffle_hesit_tile!(tiledata::HesiTileData, rng::AbstractRNG)
    shuffle!(rng,tiledata.tile_h)
end

function calc_tile_pm_hes_lin(hesitile::HesiTileData, desired_mean_pmis::AbstractFloat)
    
   AV_upm =  @.( hesitile.people * hesitile.tile_h) / sum(hesitile.people)
   Mscale = desired_mean_pmis / AV_upm
   pmis_tile_order_vhes = @. hesitile.tile_h * Mscale

   pmis_tile_order_vhes
end

function calc_pmisinf_rescale_pvhe(datatile::TileData, tiles_vhes_to_filt::AbstractDataFrame, total_people::Integer, desired_mean_pmis::AbstractFloat; col_vaxhes="pvhe")
    df = filter("tile_id" => x -> in(x, datatile.tiles_idcs), tiles_vhes_to_filt)
    sort!(df, "tile_id")
    
    x = replace(df[!,col_vaxhes], missing => NaN)
    @assert any(isnan.(x)) == false
    ## rescale
    upm = (x.-minimum(x))./(maximum(x)-minimum(x))
    ## calculate average unscaled_pm
    AV_upm = sum(datatile.tiles_pop[tid]*y for (tid,y) in zip(df[!,"tile_id"], upm) ) / total_people
    Mscale = desired_mean_pmis / AV_upm
    pmis_tile_order_vhes =  upm * Mscale
    return DataFrame(tile_id=df[!,"tile_id"], p_mis=pmis_tile_order_vhes)
end

function max_global_pmis_rescale_pvhe(datatile::TileData, tiles_vhes_to_filt::AbstractDataFrame, total_people::Integer; col_vaxhes="pvhe")
    df = filter("tile_id" => x -> in(x, datatile.tiles_idcs), tiles_vhes_to_filt)
    sort!(df, "tile_id")
    
    x = replace(df[!,col_vaxhes], missing => NaN)
    @assert any(isnan.(x)) == false
    ## rescale
    upm = (x.-minimum(x))./(maximum(x)-minimum(x))
    ## calculate average unscaled_pm
    AV_upm = sum(datatile.tiles_pop[tid]*y for (tid,y) in zip(df[!,"tile_id"], upm) ) / total_people

    return AV_upm
end