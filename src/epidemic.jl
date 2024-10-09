using ProgressMeter

using GraphEpidemics
import GraphEpidemics as GE

struct SMIRModel
    beta::AbstractFloat 
    gamma::AbstractFloat
    a::AbstractFloat ## how more probable is infection from M
    b::AbstractFloat ## how many times is more likely for M to be infected
end

function epidemic_model_SMIR(glob_model, ismis, Npop)
    betas = fill(convert(Float32,glob_model.beta), Npop)
    sigmas = ones(Float32, Npop)
    
    betas[ismis] .*= glob_model.a
    sigmas[ismis] .= glob_model.b
    
    SIRModelSus(betas, sigmas, glob_model.gamma)
end


struct SMIRStatesCounter{T<:AbstractVector} <: GraphEpidemics.AbstractStatesCounter
    ismis::T
end

import GraphEpidemics: count_states

function count_states(counter::SMIRStatesCounter,states::AbstractVector)
    stateM = states[counter.ismis]
    all = (StatsBase.counts(states,3))
    Mc = (StatsBase.counts(stateM,3))
    Oc = all .- Mc
    SVector{6, Int}(vcat(Oc,Mc))
end



function stats_misinf_tile(ismisinf::Union{Vector,BitVector}, datatile::TileData)
    [(; :tile_id => k, :n_misinf => sum(ismisinf[t]) ) for (k,t) in datatile.tiles_idcs]
end

function draw_misinformed_homogen_stats(betaprob::AbstractFloat, beta_mult_mis::Real, p_misinf::AbstractFloat,
    N::Integer,  rng::AbstractRNG, datatile::TileData)
    ismisinf, probinfs = draw_misinformed_homogen(betaprob, beta_mult_mis, p_misinf, N , rng)
    stats = stats_misinf_tile(ismisinf, datatile)

    return ismisinf, probinfs, stats
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

function stack_sims_results(r::Vector{Tuple})
    NSIMS = length(r)
    maxt = maximum(map(x->length(x), r))
    mega_arr = zeros(Int, (maxt, 3,NSIMS));
    for s=1:NSIMS
        vv = r[s]
        tt=length(r[s])
        for i in eachindex(vv[1])
            mega_arr[1:tt,i,s] = getindex.(vv,i)
            mega_arr[tt+1:end,i,s].= mega_arr[tt,i,s]
        end
        
    end
    mega_arr
end

function stack_sims_results(r::Vector{Matrix})
    NSIMS = length(r)
    maxt = maximum(map(x->size(x,1), r))
    #println(maxt)
    mega_arr = zeros(eltype(r[1]), (maxt, size(r[1],2),NSIMS));
    for s=1:NSIMS
        vv = r[s]
        tt=size(r[s],1)
        #for i in eachindex(vv[1])
        mega_arr[1:tt,:,s] = vv
        for l=1:size(vv,2)
            mega_arr[tt+1:end,l,s].= mega_arr[tt,l,s]
        end
        #end
        
    end
    mega_arr
end

function run_single_sim_city(G::AbstractGraph, T::Int, misinf_tile_dict::Dict, datatile::TileData,  epid_pars::GE.SIRModel, rng::AbstractRNG,
    kind::Symbol, num_seeds::Int, beta_mult_mis::Real, frac_misinf_equal::AbstractFloat; beta_IorS=:S)
    N = nv(G)
    #=if kind == :no_misinf 
        ismisinf, probinfs = draw_misinformed_homogen(probs_adjusted.beta, beta_mult_mis, 0.0, N , rng) ## no misinformed
        #stats = stats_misinf_tile(ismisinf, datatile)
    elseif  kind == :misinf_equal 
        ismisinf, probinfs = draw_misinformed_homogen(probs_adjusted.beta, beta_mult_mis, frac_misinf_equal, N, rng) ## avg misinformed
        #stats = stats_misinf_tile(ismisinf, datatile)
    elseif  kind == :misinf_tile
        ismisinf, probinfs, stats=draw_misinformed_tile(probs_adjusted.beta, beta_mult_mis, misinf_tile_dict, N, datatile.tiles_idcs, rng)
    else
        throw(ArgumentError("The kind argument $kind is invalid. Choose between :no_misinf, :misinf_equal, or :misinf_tile"))
    end
    =#

    
    model = GE.SIRModel(probinfs, epid_pars.gamma)

    pz = rand(rng, 1:N, num_seeds)
    #tidstart, pz = draw_people_same_tile(tids_large_pop, tiles_idcs, num_seeds, rng)
    dat, counts = GE.run_sir_fast(G, model,T, rng, pz, beta_IorS=beta_IorS)

    return counts, dat, ismisinf #DataFrame(stats)
end


function run_epidemics_city(G::AbstractGraph, T::Int, NSIMS::Int, misinf_tile_dict::Dict, datatile::TileData,  model_probs::GE.SIRModel, rng::AbstractRNG; 
                    kind::Symbol=:misinf_tile, num_seeds::Int=20)
    simdatas = GE.SIRSimData[]
    statesrec = Vector{NTuple{3,Int}}[]
    #misinfstats = DataFrame[]
    ismisinf_all = Vector{BitVector}(undef,NSIMS)
    
    #lk =ReentrantLock()
    frac_mean_mis = mean_misinf_pop(misinf_tile_dict, datatile)
    betaprob = model_probs.beta
    gprob = model_probs.gamma
    N =  nv(G)
        
    @showprogress for i=1:NSIMS
        if kind == :no_misinf 
            ismisinf, probinfs = draw_misinformed_homogen(betaprob, beta_mult_mis, 0.0, N , rng) ## no misinformed
            #stats = stats_misinf_tile(ismisinf, datatile)
        elseif  kind == :avg_misinf 
            ismisinf, probinfs = draw_misinformed_homogen(betaprob, beta_mult_mis, frac_mean_mis, N, rng) ## avg misinformed
            #stats = stats_misinf_tile(ismisinf, datatile)
        elseif  kind == :misinf_tile
            ismisinf, probinfs, stats=draw_misinformed_tile(betaprob, beta_mult_mis, misinf_tile_dict, N,datatile.tiles_idcs, rng)
        else
            throw(ArgumentError("The kind argument $kind is invalid. Choose between :no_misinf, :avg_misinf, or :misinf_tile"))
        end
        
        model = GE.SIRModel(probinfs, gprob)
    
        pz = rand(rng, 1:nv(G), num_seeds)
        #tidstart, pz = draw_people_same_tile(tids_large_pop, tiles_idcs, num_seeds, rng)
        dat, counts = GE.run_sir_fast(G, model,T, rng, pz, beta_IorS=:S)
        #lock(lk) do
            push!(statesrec, counts)
            push!(simdatas, dat)
            #print("$i  ")
            #push!(misinfstats, DataFrame(stats))
            ismisinf_all[i] = ismisinf
        #end
        #next!(p)
    end
    println("Finished")
    
    #=nstates = reshape(reduce(vcat,statesrec),T+1,NSIMS,3)#calc_nstates_all(infects,recovs, T);
    nstates = permutedims(nstates,(1,3,2))
    nstates[:,:,1];
    =#
    nstates = stack_sims_results(statesrec)

    (simdatas, nstates, ismisinf_all)
end


function run_epidemics_city_parallel(G::AbstractGraph, T::Int, NSIMS::Int, misinf_tile_dict::Dict, datatile::TileData,  model_probs::GE.SIRModel,seed::Int, a; 
                    kind::Symbol=:misinf_tile, num_seeds::Int=20, beta_mult_mis::Real=2, beta_IorS=:S)
    #simdatas = GE.SIRSimData[]
    #statesrec = Vector{NTuple{3,Int}}[]
    #misinfstats = DataFrame[]
    ismisinf_all = Vector{BitVector}(undef,NSIMS)

    
    resulock =ReentrantLock()
    frac_mean_mis = mean_misinf_pop(misinf_tile_dict, datatile)
    nthreads = Threads.nthreads()
    #RNGS = [Xoshiro(i+seed-1) for i=1:nthreads]
    #locksrng = 

    ## Use static thread assignment so that no other thread can use the same generator
    ## Slower, but safer
    #progr = Progress(NSIMS)
    simdatas = Vector{GE.SIRSimData}(undef,N)
    statesrec = Vector{Matrix}(undef, N)
    fin = 0

    Threads.@threads for i=1:NSIMS
        thid = Threads.threadid()
        #print("th $thid for $i;  ")
        rng = Xoshiro(i+seed-1)

        counts, dat, ismisinf = run_single_sim_city(G, T, misinf_tile_dict, datatile, model_probs, rng, kind, num_seeds, beta_mult_mis, frac_mean_mis; beta_IorS=beta_IorS)
        lock(resulock) do 

            statesrec[i] = counts
            simdatas[i] = dat
            fin += 1
            print("\r$fin / $NSIMS     ")

            ismisinf_all[i] = ismisinf
            #push!(misinfstats, statsdf)
            #next!(progr)
        end
    end
    #finish!(progr)
    println("\nFinished")
    
    nstates = stack_sims_results(statesrec)

    (simdatas, nstates, ismisinf_all)
end

function run_epidemics_city_SIR_parallel(G::AbstractGraph, T::Int, NSIMS::Int, datatile::TileData,  model_probs::GE.SIRModel,seed::Int; 
                    num_seeds::Int=20)
    simdatas = GE.SIRSimData[]
    statesrec = Vector{NTuple{3,Int}}[]
    ismisinf_all = Vector{BitVector}(undef,NSIMS)

    
    resulock =ReentrantLock()
    frac_mean_mis = 0.0 #mean_misinf_pop(misinf_tile_dict, datatile)
    nthreads = Threads.nthreads()

    ## Use static thread assignment so that no other thread can use the same generator
    ## Slower, but safer
    #progr = Progress(NSIMS)
    fin = 0

    Threads.@threads for i=1:NSIMS
        thid = Threads.threadid()
        #print("th $thid for $i;  ")
        rng = Xoshiro(i+seed-1)

        counts, dat, ismisinf = run_single_sim_city(G, T, Dict{String,String}(), datatile, model_probs, rng, :no_misinf, num_seeds, 0., frac_mean_mis)
        lock(resulock) do 
            push!(statesrec, counts)
            push!(simdatas, dat)
            fin += 1
            print("\r$fin / $NSIMS     ")
            #next!(progr)
            ismisinf_all[i] = ismisinf
        end
    end
    #finish!(progr)
    println("\nFinished")
    
    nstates = stack_sims_results(statesrec)

    (simdatas, nstates, ismisinf_all)
end

