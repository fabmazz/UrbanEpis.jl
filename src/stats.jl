vec_mean_std(x) = [mean_and_std(x)]
empty_bitvector(N::Integer) = BitVector(zeros(Bool, N))

function calc_num_connection_tiles(G::SimpleGraph, tile_for_i::AbstractVector)
    tiles_conns = Dict{Tuple{Int,Int}, Int}()
    for e in (edges(G))
        t_i = tile_for_i[e.src]
        t_j = tile_for_i[e.dst]
        for c in [ (t_i,t_j), (t_j,t_i)]
            if !(c in keys(tiles_conns))
                tiles_conns[c] = 1
            else
                tiles_conns[c]+=1
            end
        end
    end
    tiles_conns
end

function convert_dict_to_df(mdict::Dict, keys_names, val_name::String)
    df=DataFrame(Tables.table(stack(keys(mdict))'),keys_names)
    df[!,val_name] = collect(values(mdict))
    df
end

function calc_flow_infection_tiles(simd::GraphEpidemics.SIRSimData, tiles_conns::Dict, tile_for_i::AbstractVector)
    du_out = Dict(k=> 0 for k in keys(tiles_conns));
    N = size(simd.infect_node,1)

    for i=1:N
        if isnan(simd.infect_time[i])
            continue
        end
        tid = tile_for_i[i]
        if simd.infect_node[i] > 0
            infect_tid = tile_for_i[Int(simd.infect_node[i])]
            du_out[(infect_tid, tid)] += 1
        end
    end

    df=convert_dict_to_df(du_out, ["tile_from","tile_to"], "num_infect")
    sort!(df, ["tile_from","tile_to"])
end

function join_all_flow_infect(dataframes::Vector{<:DataFrame})
    if length(dataframes)<2
        return dataframes[1]
    end
    
    refdf = dataframes[1]
    sort!(refdf,["tile_from","tile_to"] )
    tif = refdf[!,"tile_from"]
    tito = refdf[!,"tile_to"]
    dfout = DataFrame(Dict("tile_from"=>tif, "tile_to"=>tito, "ninf_1" => refdf[!,"num_infect"]))
    for i=2:length(dataframes)
        d = dataframes[i]
        sort!(d, ["tile_from","tile_to"])
        @assert all(d[!,"tile_from"].==tif)
        @assert all(d[!,"tile_to"] .== tito)
        dfout[!,"ninf_$i"] = d[!,"num_infect"]
    end

    dfout
end

function save_trace_inf_misinf(simdata::GraphEpidemics.SIRSimData, ismis, counts)
    ##counts is a matrix of shape [T,N_states=6]
    N = length(simdata.infect_node)
    #st = zeros(Int, ())
    T = size(counts,1)
    for i = 1:N
        if ismis[i]
            offs=3
        else
            offs=0
        end

        if isnan(simdata.infect_time[i]) || isinf(simdata.infect_time[i])
            ## never infected
            counts[1:T, 1+offs] .+= 1
        else
            tinf = convert(Int, simdata.infect_time[i])+1
            ## offset, first time is actually 0 -> 0:tinf-1 => 1:tinf
            tinf_act = max(tinf+1,1)
            counts[1:(tinf_act-1), 1+offs] .+=1
            # 0:trec  -> 1:trec+1
            trec = tinf+simdata.rec_delays[i] +1
            #trec_act = Lmin(trec,)
            if trec > T
                counts[tinf_act:T, 2+offs].+=1
            else
                counts[tinf_act:trec-1,2+offs].+=1
                counts[trec:T, 3+offs].+=1
            end
        end
    end
end
function get_times_infs_recs(simd::GE.SIRSimData, idcs, tstep::AbstractFloat=1.0)
    l=@. ~isnan(simd.infect_time[idcs])
    iinf = idcs[l]
    iis= @. (simd.infect_time[iinf] +1)*tstep
    dels = (simd.rec_delays[iinf] .*tstep) .+ iis
    iinf,iis, dels
end

function calc_tiles_inf_history(datatile::TileData, simd::GE.SIRSimData, misbeh::AbstractVector, nTimesteps::Integer)
    fracinf = zeros(nTimesteps, length(unique(datatile.tiles_idcs)),3)
    for (it,tileid) in enumerate(datatile.tiles_idcs)
        iinf, ti,tr=get_times_infs_recs(simd, datatile.idcs_in_tile[tileid], 1.)
        
        frinf = length(iinf) / length(datatile.idcs_in_tile[tileid])
        
        Ntile = StatsBase.counts(misbeh[datatile.idcs_in_tile[tileid]])
        ismis_tile = misbeh[iinf];
        
        #ds = DataFrame(:idcs => iinf, :tinf => ti, :trec => tr, :ismis => ismis_tile)
        cstats=zeros(nTimesteps,2);
        
        for i=eachindex(tr)
            tr_sel = min(tr[i]+1, size(cstats,1))
            u = ismis_tile[i] ? 2 : 1
            cstats[convert(Int,ti[i]+1):convert(Int,tr_sel), u] .+= 1
        end
        
        fracinf[:,it,2:3] = cstats./Ntile'
        fracinf[:,it,1] = sum(cstats, dims=2)[:,1] / sum(Ntile)
    end
    fracinf
end

function find_peak_t_values_fracinf(fracinf::AbstractArray, tstep::Real)
    vs, idx = findmax(fracinf, dims=1)

    @assert all(map(u-> u[2], idx[1,:,1]) == collect(1:size(fracinf,2)) )
    times_pk_all=stack(map(u-> (u[1]-1)*tstep, idx[1,:,k]) for k=1:3)
    peaks_i_all = vs[1,:,:]

    times_pk_all, peaks_i_all
end

"""
`calc_tiles_timeinf_infector(datatile::TileData,simd::GE.SIRSimData, tiles_for_i::AbstractVector)`

Calculate the time of first infection and the time of last recovery in a tile, and the tile index from which the infection came from.
Return a DataFrame with the columns "tile_id", "time_inf", "tile_inf", "tile_rec".
"""
function calc_tiles_timeinf_infector(datatile::TileData,simd::GE.SIRSimData, tiles_for_i::AbstractVector)
    N = size(simd.infect_node,1)
    u=DataFrame(:id=>1:N,:tile_id => tiles_for_i[1:N], :time_inf => simd.infect_time .+1, :infector => simd.infect_node)
    u[!,:time_rec] = simd.rec_delays .+ u.time_inf
    u = filter(:time_inf => x-> !isnan(x), u )
    u[!,:tile_inf] = map(i-> i>0 ? Int(tiles_for_i[Int(i)]) : -10, u.infector)
    
    maxrec = combine(groupby(u, :tile_id), :time_rec => maximum => :tile_rec)
    maxrec[!,:tile_rec] = convert.(Int, maxrec.tile_rec)

    uextra=u[u.tile_id .!= u.tile_inf,:]
    result = combine(groupby(uextra, :tile_id)) do group
        min_row = argmin(group.time_inf) # Find the index of the minimum y in each group
        (time_inf = Int(group.time_inf[min_row]), tile_inf=group.tile_inf[min_row])
    end
    result = innerjoin(result, maxrec, on="tile_id")
    #z=fill(-10,(Ntiles, 2))
    iis = setdiff(datatile.tiles_idcs,result.tile_id)
    resu=vcat(result,DataFrame(tile_id = iis, time_inf=-50, tile_inf=-100, tile_rec=-50))
    sort!(resu, :tile_id)

    @assert resu.tile_id == datatile.tiles_idcs

    resu
end


"""
`calc_frac_infected_tile_SMIR(simdat::GE.SIRSimData,datatile::TileData, ismis::Union{BitVector, Vector{Bool}})`

Calculate the fraction of people (both O and M) that are infected in each tile.

Takes the `SIRSimData` object from the simulation, the `TileData` object which contains the information on the tiles,
and the vector indicating which individuals are Misbehaving

"""
function calc_frac_infected_tile_SMIR(simdat::GE.SIRSimData,datatile::TileData, ismis::Union{BitVector, Vector{Bool}})
    ntiles = size(datatile.tiles_idcs,1)
    infect_tiles = Matrix{Float32}(undef,ntiles,3)
    NM_t = Vector{Int}(undef,ntiles)
    # = zeros(size(datatile.tiles_idcs))
    inf_all = @. !isnan(simdat.infect_time)
    for i in eachindex(datatile.tiles_idcs)
        tid = datatile.tiles_idcs[i]
        peop = datatile.idcs_in_tile[tid]
        Np = size(peop,1)
        infect = @view inf_all[peop]
        mi = @view ismis[peop]
        nM = sum(mi)
        mandinf = @view infect[mi]
        oandinf = @.(infect & !mi) #@view infect[.!mi]
        infect_tiles[i,2] = sum(mandinf) / (nM)
        infect_tiles[i,1] = sum(oandinf) / (Np-nM)
        infect_tiles[i,3] = sum(infect) / Np
        NM_t[i] = nM
    end

    infect_tiles, NM_t
end

function count_SMIR_state(states)
    c = SVector{6,Int}(
        sum(states.==i) for i=1:6
    )
end

function calc_tile_infect_stats(simdat, isM, datatile)
    resu = NamedTuple[]
    for (tid, ids) in datatile.tiles_idcs
        misinfs = isM[ids]
        times_ids = simdat.infect_time[ids]
        infs_ids = @. ~isnan(times_ids)
        nM  = sum(misinfs)
        
        #println(misinfs," ", infs_ids)
        totinf_tile = sum(infs_ids)
        ninfs_M_tile = sum(infs_ids[misinfs])
        ar_all = totinf_tile / length(ids)
        ar_m = nM > 0 ? ninfs_M_tile / nM : 0
        ar_notm = (totinf_tile - ninfs_M_tile) / (length(ids) - nM)
        #println("totinf $totinf_tile -> $(round_r(ar_all,3)),\t I_M = $ninfs_M_tile -> $(round_r(ar_m,3))\t I_0 = $(totinf_tile - ninfs_M_tile) -> $(round_r(ar_notm,3)) ")
    
        push!(resu, (; :tile_id => tid, :frac_inf => ar_all, :ninf_m => ninfs_M_tile, :n_misinf => nM, :ninf_o => (totinf_tile - ninfs_M_tile), :n_o => (length(ids) - nM)))
    end
    
    df_res = DataFrame(resu)
    
    ar_global = mean( @. ~isnan(simdat.infect_time))
    df_res[!,"rel_att"] = df_res[!,:frac_inf] / ar_global

    df_res
end

covar(a,b) = mean((a.-mean(a)) .* (b.-mean(b)))

function find_avg_attack_rates_allsims(simdatas, misinfsarr, datatile)
    NSIMS = length(simdatas)

    all_fr=[calc_tile_infect_stats(simdatas[m], misinfsarr[m], datatile) for m=1:NSIMS]

    all_fr_cat = reduce(vcat, all_fr)

    cols_m = ["ninf_m", "n_misinf", "ninf_o", "n_o"]

    ops=[c => vec_mean_std => ["$(c)_mean", "$(c)_std"] for c in cols_m]
    covars = [["n_misinf", "ninf_m"] => covar => "M_inf_n_covar", ["ninf_o", "n_o"] => covar => "O_inf_n_covar"]
    combine(groupby(all_fr_cat,:tile_id), "frac_inf" => vec_mean_std => ["mean_finf", "finf_std"], 
            "rel_att" => vec_mean_std => ["rel_att_mean", "rel_att_std"],ops..., covars...
        )
end

### DEPRECATED
#=
function find_avg_attack_rates_alldata_legacy(simdatas, datatile)
    alldas = DataFrame[] #s=simdatas[1]
    minfs=[]
    for s in simdatas
        dd=find_fraction_infected((@. ~isnan(s.infect_time)), datatile)
        infected = @. ~isnan(s.infect_node)
        minf = mean(infected)
        frac_inf = collect(values(dd))
        rel_ar = frac_inf ./ minf
        push!(minfs, minf)
        df=DataFrame(:tile_id => collect(keys(dd)), :frac_inf => frac_inf, :rel_att => rel_ar)
        push!(alldas,df)
    end
    dff = reduce(vcat,alldas)
    combine(groupby(dff,:tile_id), "frac_inf" => vec_mean_std => ["mean_finf", "finf_std"], "rel_att" => vec_mean_std => ["rel_att_mean", "rel_att_std"]), minfs
end

function tiles_misinfect_stats_full(misinfstats, tile_data; people_col=:people)
    mis_stats_shuff = reduce(vcat,misinfstats);
    mean_misinf_shuff = combine(groupby(mis_stats_shuff,:tile_id), :n_misinf => vec_mean_std => [:mean_n_misinf, :n_misinf_std], )
    mean_misinf_shuff = innerjoin(mean_misinf_shuff, tile_data, on="tile_id")
    mean_misinf_shuff[!,:frac_misinf_mean] = @. mean_misinf_shuff.mean_n_misinf / mean_misinf_shuff[!,people_col]

    mean_misinf_shuff
end
=#
function moving_average(arr::Vector, win_size::Integer)
    n = length(arr)
    avg = zeros(n)
    for i=1:n
        ms = max(1,i-win_size+1)
        avg[i] = mean(arr[ms:i])
    end
    avg
end
#=
function calc_curves_misinf_slower(simds,counts_tot, misinfidcs)
    TT=size(counts_tot,1)
    tmax=TT-1
    counts_misinf = zeros(Int,size(counts_tot));
    tarrs = collect(0:tmax)
    i=1
    for (sd, ism) in zip(simds, misinfidcs)

        dels = sd.rec_delays
        infts = sd.infect_time.+1
        trecs = infts .+ dels

        trM = trecs[ism]
        tiM = infts[ism]
        
        isInf =(@. tarrs >= tiM') 
        counts_misinf[:,1,i] .= vec(sum(.!isInf, dims=2))
        nInfected=vec(sum(isInf, dims=2))
        counts_misinf[:,3,i] .= vec(sum((@. tarrs >= trM'),dims=2))
        counts_misinf[:,2,i] .= nInfected .-counts_misinf[:,3,i]

        i+=1
    end

    counts_misinf
end

function calc_curves_misinf(simds,counts_tot, misinfidcs)
    TT=size(counts_tot,1)
    tmax=TT-1
    counts_misinf = zeros(Int,size(counts_tot));
    tarrs = collect(0:tmax)
    NSIMS = length(simds)
    for i=1:NSIMS
        sd = simds[i]
        ism = misinfidcs[i]

        dels = sd.rec_delays
        infts = sd.infect_time.+1
        trecs = infts .+ dels

        trM = trecs[ism]
        tiM = infts[ism]
        num_M = sum(ism)
        
        for t in tarrs
            isinf=@. t >= tiM
            ninfect = sum(isinf)
            nR = sum(@. t >= trM)
            nI = ninfect - nR
            nS = sum(.!isinf) #num_M - ninfs
            ti=t+1
            counts_misinf[ti,1,i] = nS
            counts_misinf[ti,2,i] = nI
            counts_misinf[ti,3,i] = nR
        end

    end

    counts_misinf
end

function calc_curves_misinf_parallel(simds,counts_tot, misinfidcs)
    TT=size(counts_tot,1)
    tmax=TT-1
    counts_misinf = zeros(Int,size(counts_tot));
    tarrs = collect(0:tmax)
    NSIMS = length(simds)
    
    @Threads.threads for i=1:NSIMS
        sd = simds[i]
        ism = misinfidcs[i]

        dels = sd.rec_delays
        infts = sd.infect_time.+1
        trecs = infts .+ dels

        trM = trecs[ism]
        tiM = infts[ism]
        #num_M = sum(ism)
        
        for t in tarrs
            isinf=@. t >= tiM
            ninfect = sum(isinf)
            nR = sum(@. t >= trM)
            nI = ninfect - nR
            nS = sum(.!isinf) #num_M - ninfs
            ti=t+1
            counts_misinf[ti,1,i] = nS
            counts_misinf[ti,2,i] = nI
            counts_misinf[ti,3,i] = nR
        end

    end

    counts_misinf
end
=#
