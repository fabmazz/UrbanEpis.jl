vec_mean_std(x) = [mean_and_std(x)]

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
    combine(groupby(dff,:tile_id), :frac_inf => vec_mean_std => [:mean_finf, :finf_std], :rel_att => vec_mean_std => [:rel_att_mean, :rel_att_std]), minfs
end

function tiles_misinfect_stats_full(misinfstats, tile_data; people_col=:people)
    mis_stats_shuff = reduce(vcat,misinfstats);
    mean_misinf_shuff = combine(groupby(mis_stats_shuff,:tile_id), :n_misinf => vec_mean_std => [:mean_n_misinf, :n_misinf_std], )
    mean_misinf_shuff = innerjoin(mean_misinf_shuff, tile_data, on="tile_id")
    mean_misinf_shuff[!,:frac_misinf_mean] = @. mean_misinf_shuff.mean_n_misinf / mean_misinf_shuff[!,people_col]

    mean_misinf_shuff
end