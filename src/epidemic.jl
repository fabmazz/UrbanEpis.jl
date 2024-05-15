using ProgressMeter

import GraphEpidemics as GE

function stats_misinf_tile(ismisinf::Union{Vector,BitVector}, datatile::TileData)
    [(; :tile_id => k, :n_misinf => sum(ismisinf[t]) ) for (k,t) in datatile.tiles_idcs]
end

function mean_misinf_pop(misinf_tile_dict::Dict, datatile::TileData)
    c = 0
    s = 0.0
    for (k,fr) in misinf_tile_dict
        s += datatile.tiles_pop[k] * fr
        c += datatile.tiles_pop[k]
    end
    s/c
end

function run_epidemics_city(G::AbstractGraph, T::Int, NSIMS::Int, misinf_tile_dict::Dict, datatile::TileData,  model_probs::GE.SIRModel, rng::AbstractRNG; 
                    kind::Symbol=:misinf_tile, num_seeds::Int=20, beta_mult_mis::Real=2)
    simdatas = GE.SIRSimData[]
    statesrec = Matrix{Int}[]
    misinfstats = DataFrame[]
    
    #lk =ReentrantLock()
    frac_mean_mis = mean_misinf_pop(misinf_tile_dict, datatile)
    betaprob = model_probs.beta
    gprob = model_probs.gamma
    N =  nv(G)
        
    @showprogress for i=1:NSIMS
        if kind == :no_misinf 
            ismisinf, probinfs = draw_misinformed_homogen(betaprob, beta_mult_mis, 0.0, N , rng) ## no misinformed
            stats = stats_misinf_tile(ismisinf, datatile)
        elseif  kind == :avg_misinf 
            ismisinf, probinfs = draw_misinformed_homogen(betaprob, beta_mult_mis, frac_mean_mis, N, rng) ## avg misinformed
            stats = stats_misinf_tile(ismisinf, datatile)
        elseif  kind == :misinf_tile
            ismisinf, probinfs, stats=draw_misinformed_tile(betaprob, beta_mult_mis, misinf_tile_dict, N,datatile.tiles_idcs, rng)
        else
            throw(ArgumentError("The kind argument $kind is invalid. Choose between :no_misinf, :avg_misinf, or :misinf_tile"))
        end
        
        model = GE.SIRModel(probinfs, gprob)
    
        pz = rand(rng, 1:nv(G), num_seeds)
        #tidstart, pz = draw_people_same_tile(tids_large_pop, tiles_idcs, num_seeds, rng)
        dat, counts = GE.run_sir_fast(G, model,T, rng, pz,prob_infect_I=false)
        #lock(lk) do
            push!(statesrec, counts)
            push!(simdatas, dat)
            #print("$i  ")
            push!(misinfstats, DataFrame(stats))
        #end
        #next!(p)
    end
    println("Finished")
    
    nstates = reshape(reduce(vcat,statesrec),T+1,NSIMS,3)#calc_nstates_all(infects,recovs, T);
    nstates = permutedims(nstates,(1,3,2))
    nstates[:,:,1];

    (simdatas, nstates, misinfstats)
end