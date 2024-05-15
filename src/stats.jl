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