vec_mean_std(x) = [mean_and_std(x)]

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

function moving_average(arr::Vector, win_size::Integer)
    n = length(arr)
    avg = zeros(n)
    for i=1:n
        ms = max(1,i-win_size+1)
        avg[i] = mean(arr[ms:i])
    end
    avg
end