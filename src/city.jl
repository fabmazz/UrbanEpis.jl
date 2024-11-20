struct TileData{I<:Integer}
    idcs_in_tile::Dict{I,Vector{I}}
    tiles_idcs::Vector{I}
    tiles_pop::Dict{I,I}
    tile_for_i::Dict{I,I}
end

function extract_tile_data(nodes_dat_new::DataFrame; col_i="i")

    idcs_in_tile = Dict{Int,Vector{Int}}()
    tiles_pop = Dict{Int,Int}()
    tile_for_i = Dict{Int,Int}()
    for r in eachrow(nodes_dat_new)
        if r[col_i] == -1
            continue
        end
        inew = r[col_i]
        tid = r["tile_id"]
        if !(tid in keys(idcs_in_tile))
            idcs_in_tile[tid] = [inew]
            tiles_pop[tid] = 1
        else
            push!(idcs_in_tile[tid],inew)
            tiles_pop[tid]+=1
        end
        tile_for_i[inew] = tid
    end

    TileData(idcs_in_tile,sort(collect(keys(idcs_in_tile))), tiles_pop, tile_for_i)
end

tile_for_pop(tiledata::TileData, i::Int) = tiledata.tile_for_i[i]

function get_tiles_with_pop(data::TileData,condition::Function)
    idcs = Int[]
    for (k,v) in data.tiles_pop
        if condition(v)
            push!(idcs, k)
        end
    end
    idcs
end

function draw_people_same_tile(sample_tiles_idcs::Vector, tiles_idcs::Dict{I,Vector{I}}, npeople::Integer, rng::AbstractRNG) where I <: Integer
    tid = rand(rng, sample_tiles_idcs)
    idcs = rand(rng, tiles_idcs[tid], npeople)
    tid,idcs
end

function draw_misinformed(betan::AbstractFloat, multi_misinf::Real, p_misinf::AbstractFloat,
    gr::AbstractGraph,
    rng::AbstractRNG)
        ismisinf = rand(rng,nv(gr)) .< p_misinf;
        p=convert(Float32, betan)
        betasall = fill(p, nv(gr))
        betasall[ismisinf] .= p*multi_misinf

        ismisinf, betasall
end

function find_fraction_infected(isinfected::BitVector, tilesdata::TileData)
    find_fraction_infected(isinfected, tilesdata.idcs_in_tile)
end

function find_fraction_infected(isinfected::BitVector, tiles_idcs::Dict{I,Vector{I}}) where I <: Integer
    dats = Dict{I,Float64}()
    for (tid, ids) in tiles_idcs
        s = isinfected[ids]
        dats[tid] = sum(s) / length(s)
    end
    dats
end


function find_avg_frac_inf_alldata(simdatas, datatile)
    alldas = DataFrame[] #s=simdatas[1]
    for s in simdatas
        dd=find_fraction_infected((@. ~isnan(s.infect_time)), datatile)
        df=DataFrame(:tile_id => collect(keys(dd)), :frac_inf => collect(values(dd)))
        push!(alldas,df)
    end
    dff = reduce(vcat,alldas)
    combine(groupby(dff,:tile_id), :frac_inf => mean => :mean_inf)
end

convert_arr_tile_id_dict(arr::Vector, datatile::TileData) = Dict(zip(datatile.tiles_idcs, arr))

#=function draw_misinformed_tile(betan::AbstractFloat, mult::Real, p_misinf_tile::Dict, N::Integer, tile_people::Dict, rng::AbstractRNG)
    #convert(Float32, betan)
    betasall = fill(betan, N)
    ismisinf = fill(false, N)

    stats = NamedTuple[]
    for (tid, people) in tile_people
        misinf  = rand(rng, length(people)) .< p_misinf_tile[tid]
        ii = people[misinf]
        betasall[ii] .= mult*betan
        ismisinf[ii] .= true
        
        push!(stats, (; :tile_id => tid, :n_misinf => sum(misinf)))
    end

    ismisinf, betasall, stats
end
=#

function draw_misbeh_hom(betan::AbstractFloat, mult::Real, p_misinf::AbstractFloat,
            N::Integer,  rng::AbstractRNG)

    ismisinf = rand(rng,N) .< p_misinf;
    #p=convert(beta_type, betan)
    betasall = fill(betan, N)
    betasall[ismisinf] .= betan*mult

    ismisinf, betasall
end
function calc_tiles_infect_stats(tiles_idcs::Dict, infect_time::Vector{<:AbstractFloat})
    ## assuming that the indcs for each tile are the same as the graph (1:N)
    stats = NamedTuple[]
    for (tid) in sort(collect(keys(tiles_idcs)))
        idcs = tiles_idcs[tid]
        tinfects= infect_time[idcs]
        ninfects = sum(@. ~isnan(tinfects))
        frac_infect = ninfects / length(idcs)
        t = nanminimum(tinfects)+1
        push!(stats, (; :tile_id => tid, :time_infect=>convert(Float32,t), :n_inf => ninfects, :frac_inf=>frac_infect))
    end
    stats
end

function draw_misbehav_tiles(tiles_idcs::AbstractVector, pmisinf_tile::Vector, idcs_in_tile::Dict, N::Integer, rng::AbstractRNG)
    ismisinf = empty_bitvector(N)#fill(false, N)
    frac_misinf = zeros(length(tiles_idcs))
    for i=eachindex(tiles_idcs) #(t, pm) in zip(datatile.tiles_idcs, pmisinf_tile)
        tid = tiles_idcs[i]
        people = idcs_in_tile[tid]
        misinf  = rand(rng, length(people)) .< pmisinf_tile[i]
        ii = people[misinf]
        ismisinf[ii] .= true
        frac_misinf[i] = sum(misinf)/length(people)
    end
    ismisinf, frac_misinf
end