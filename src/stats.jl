function find_avg_attack_rates_alldata(simdatas, datatile)
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
    combine(groupby(dff,:tile_id), :frac_inf => mean => :mean_inf, :rel_att => mean), minfs
end

function tiles_misinfect_stats_full(misinfstats, tile_data)
    mis_stats_shuff = reduce(vcat,misinfstats);
    mean_misinf_shuff = combine(groupby(mis_stats_shuff,:tile_id), :n_misinf => mean => :mean_n_misinf)
    mean_misinf_shuff = innerjoin(mean_misinf_shuff, tile_data, on="tile_id")
    mean_misinf_shuff[!,:frac_misinf_mean] = @. mean_misinf_shuff.mean_n_misinf / mean_misinf_shuff.count

    mean_misinf_shuff
end