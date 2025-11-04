function count_epistates_by_tileid(data::GraphEpidemics.SEIRSimData, citydata::TileData, t::Real; dtype::DataType=Int)
    # Unique sorted keys
    unique_keys = citydata.tiles_idcs
    nkeys = length(unique_keys)
    
    # Define state order and mapping
    state_syms = (:S, :E, :I, :R)
    state_to_col = Dict(:S => 1, :E => 2, :I => 3, :R => 4)
    
    # Initialize count matrix (rows = keys, cols = states)
    counts = zeros(dtype, nkeys, length(state_syms))
    
    # Loop once through all individuals, compute + accumulate
    for (row,tid) in enumerate(unique_keys)
        for i in citydata.idcs_in_tile[tid]
            ##k = citydata.tile_for_i[i]
            ##row = key_to_row[k]
            
            inf_t = data.infect_time[i]
            lat = data.lat_delays[i]
            rec = data.rec_delays[i]
            
            # determine state directly
            state_col = if isnan(inf_t)
                1
            elseif t < inf_t
                1  # :S
            elseif t < inf_t + 1 + lat
                2  # :E
            elseif t < inf_t + 1 + lat + rec
                3  # :I
            else
                4  # :R
            end
            
            counts[row, state_col] += 1
        end
    end
    
    return counts
end 

not_and_y(x,y)= (!x) && y

function assign_time!(inf_t, lat, rec, t, idxt, res, j)
    state_col = if isnan(inf_t)
        1
    elseif t < inf_t
        1  # :S
    elseif t < inf_t + 1 + lat
        2  # :E
    elseif t < inf_t + 1 + lat + rec
        3  # :I
    else
        4  # :R
    end

    res[idxt,state_col,j] +=1
end



function count_states_by_tileid_times(data::GraphEpidemics.SEIRSimData, citydata::TileData, times::AbstractVector{<:Real}; dtype::DataType=Int)
    # Unique sorted keys
    unique_keys = citydata.tiles_idcs
    nkeys = length(unique_keys)
    
    # Define state order and mapping
    state_syms = (:S, :E, :I, :R)
    state_to_col = Dict(:S => 1, :E => 2, :I => 3, :R => 4)
    
    # Initialize count matrix (rows = keys, cols = states)
    ntimes = length(times)
    counts = zeros(dtype, ntimes, nkeys, length(state_syms))
    
    # Loop once through all individuals, compute + accumulate
    for (row,tid) in enumerate(unique_keys)
        for i in citydata.idcs_in_tile[tid]
            ##k = citydata.tile_for_i[i]
            ##row = key_to_row[k]
            
            inf_t = data.infect_time[i]
            lat = data.lat_delays[i]
            rec = data.rec_delays[i]
            
            # determine state directly

            #fill!(counts[])
            if isnan(inf_t)
                ### always S
                @view(counts[:,row,1]) .+= 1
            else
                mm = times .< inf_t
                @view(counts[mm, row, 1]) .+=1
                g = times .< (inf_t +1+lat)
                mm = not_and_y.(mm, g) #@. !mm && g
                @views(counts[mm, row, 2]) .+=1
                g2 = times .< (inf_t + 1+ lat+rec)
                mm = not_and_y.(g, g2)  #@. !g && g2
                @view(counts[mm, row, 3]) .+=1
                ### for R 
                @view(counts[.!g2, row, 4]) .+=1
            end
            #=state_col = if isnan(inf_t)
                1
            elseif t < inf_t
                1  # :S
            elseif t < inf_t + 1 + lat
                2  # :E
            elseif t < inf_t + 1 + lat + rec
                3  # :I
            else
                4  # :R
            end
            
            counts[row, state_col] += 1
            =#
        end
    end
    
    return counts
end 



function count_states_by_group(data::GraphEpidemics.SEIRSimData, groupForI::AbstractVector{<:Integer}, times::AbstractVector{<:Real};
     dtype::DataType=Int, fractional::Bool=false)
    # Unique sorted keys
    unique_keys = unique(groupForI)#citydata.tiles_idcs
    nkeys = length(unique_keys)
    N = length(groupForI)
    @assert length(data.infect_time) == N
    
    # Define state order and mapping
    state_syms = (:S, :E, :I, :R)
    #state_to_col = Dict(:S => 1, :E => 2, :I => 3, :R => 4)
    row_for_c = Dict(g => i for (i, g) in enumerate(unique_keys))
    
    # Initialize count matrix (rows = keys, cols = states)
    if fractional
        if dtype==Int
            dtype = Float64
        end
        @assert dtype <: AbstractFloat "Chosen dtype is not compatible for fractional (must be <: AbstractFloat)"
    end
    ntimes = length(times)
    counts = zeros(dtype, ntimes, nkeys, length(state_syms))
    
    #for (row,tid) in enumerate(unique_keys)
        for i=1:N
            row=row_for_c[groupForI[i]]
            inf_t = data.infect_time[i]
            lat = data.lat_delays[i]
            rec = data.rec_delays[i]
            
            # determine state directly

            #fill!(counts[])
            if isnan(inf_t)
                ### always S
                @view(counts[:,row,1]) .+= 1
            else
                mm = times .< inf_t
                @view(counts[mm, row, 1]) .+=1
                g = times .< (inf_t +1+lat)
                mm = not_and_y.(mm, g) #@. !mm && g
                @views(counts[mm, row, 2]) .+=1
                g2 = times .< (inf_t + 1+ lat+rec)
                mm = not_and_y.(g, g2)  #@. !g && g2
                @view(counts[mm, row, 3]) .+=1
                ### for R 
                @view(counts[.!g2, row, 4]) .+=1
            end

        end
    if fractional
        for (i,g) in enumerate(unique_keys)
            nn = sum(groupForI .== g)
            counts[:,i,:] ./= nn
        end
    end
    return counts
end 

    


function count_states_by_tileid_times_f(data::GraphEpidemics.SEIRSimData, citydata::TileData, times::AbstractVector{<:Real}; dtype::DataType=Int)
    # Unique sorted keys
    unique_keys = citydata.tiles_idcs
    nkeys = length(unique_keys)
    
    # Define state order and mapping
    state_syms = (:S, :E, :I, :R)
    state_to_col = Dict(:S => 1, :E => 2, :I => 3, :R => 4)
    
    # Initialize count matrix (rows = keys, cols = states)
    ntimes = length(times)
    counts = zeros(dtype, ntimes, length(state_syms), nkeys)
    itimes = collect(1:ntimes)
    
    # Loop once through all individuals, compute + accumulate
    for (row,tid) in enumerate(unique_keys)
        for i in citydata.idcs_in_tile[tid]
            ##k = citydata.tile_for_i[i]
            ##row = key_to_row[k]
            
            inf_t = data.infect_time[i]
            lat = data.lat_delays[i]
            rec = data.rec_delays[i]
            
            assign_time!.(inf_t, lat, rec, times,itimes,counts, row)
        end
    end
    
    return counts
end 



function count_bins!(vals::AbstractVector, bin_edges::Union{AbstractVector, AbstractRange}, bin_adds::AbstractVector)
    vsort = sort(vals)
    N = length(vals)
    c=1#vsort[1]
    idx_startp = findfirst(bin_edges .> vsort[1])
    if isnothing(idx_startp)
        ### do not have to add anything
        return
    end
    #while (idx_startp )
    if idx_startp > 1
        idx_startp -=1
    end
    idxbins=idx_startp
    #println("First index is $idxbins")
    #i_arr_s  = findfirst(vsort .>= bin_edges[idx_startp])
    for i=1:N
        val = (vsort[i])
        #print(val, "  ", bin_edges[idxbins], "  ", bin_edges[idxbins+1])
        if val < bin_edges[idxbins]
            #println(" C")
            continue
        elseif val >= last(bin_edges)
            ### we are out of the bins
            #println(" Break")
            break
        end
        while (val >= bin_edges[idxbins+1])
            ##advance
            idxbins+=1
            #println(" Advance")
            if val >= last(bin_edges)
                break
            end
        end
        
        bin_adds[idxbins] +=1
        #println(" ")
    end
    bin_adds
end

function count_bins(vals::AbstractVector, bin_edges::Union{AbstractVector, AbstractRange})
    bins_v = zeros(Int,length(bin_edges)-1)
    count_bins!(vals, bin_edges, bins_v)
end


function bin_counts_tileids(A::Vector{Vector{T}}, bin_edges::AbstractVector{<:AbstractFloat}, tileData::TileData) where T
    # Il numero di righe nell'array A
    n_rows = length(A)
    # Il numero di bin sarà il numero di bordi - 1
    n_bins = length(bin_edges) - 1
    nkeys = length(tileData.tiles_idcs)
    # Pre-alloca la matrice dei risultati per le performance
    # Risultato è una matrice n_rows x n_bins di Int
    counts_matrix = zeros(Int, n_rows, nkeys,n_bins)
    
    # 1. Definisce l'oggetto istogramma con i bordi desiderati
    # closed=:left significa che i bin sono intervalli chiusi a sinistra: [a, b)
    #h = Histogram(bin_edges, closed=:left)
    
    # 2. Itera su ogni riga dell'array 2D
    for (t, row) in enumerate(A)
        # fit(h, row_view) calcola l'istogramma solo sui dati della riga corrente.
        # .weights contiene i conteggi per ogni bin.
        for k =1:nkeys
            tid = tileData.tiles_idcs[k]
            idcs_in = tileData.idcs_in_tile[tid]

            #counts_matrix[t, k,:] = fit(Histogram,@view(row[idcs_in]), bin_edges, closed=:left).weights
            count_bins!(@view(row[idcs_in]), bin_edges,@view(counts_matrix[t,k,:]))
        end
    end
    
    return counts_matrix
end

_mean(x::AbstractVector) = sum(x)/length(x)

function mean_counts_tileids(A::Vector{Vector{T}}, tileData::TileData) where T
    # Il numero di righe nell'array A
    n_rows = length(A)
    # Il numero di bin sarà il numero di bordi - 1
    #n_bins = length(bin_edges) - 1
    nkeys = length(tileData.tiles_idcs)
    # Pre-alloca la matrice dei risultati per le performance
    # Risultato è una matrice n_rows x n_bins di Int
    means_matr = zeros(n_rows, nkeys)
    
    # 1. Definisce l'oggetto istogramma con i bordi desiderati
    # closed=:left significa che i bin sono intervalli chiusi a sinistra: [a, b)
    #h = Histogram(bin_edges, closed=:left)
    
    # 2. Itera su ogni riga dell'array 2D
    for (t, row) in enumerate(A)
        # .weights contiene i conteggi per ogni bin.
        for k =1:nkeys
            tid = tileData.tiles_idcs[k]
            idcs_in = tileData.idcs_in_tile[tid]

            #counts_matrix[t, k,:] = fit(Histogram,@view(row[idcs_in]), bin_edges, closed=:left).weights
           # count_bins!(@view(row[idcs_in]), bin_edges,@view(means_matr[t,k,:]))

            means_matr[t,k] = _mean(@view(row[idcs_in]))
        end
    end
    
    return means_matr
end

function count_indiv_prop_per_tile(indiv_prop::AbstractVector, tileData::TileData, fun::Function; dtype=Float64)
    nkeys = length(tileData.tiles_idcs)
    res = zeros(dtype,nkeys)
        for k =1:nkeys
        tid = tileData.tiles_idcs[k]
        idcs_in = tileData.idcs_in_tile[tid]

        #counts_matrix[t, k,:] = fit(Histogram,@view(row[idcs_in]), bin_edges, closed=:left).weights
        # count_bins!(@view(row[idcs_in]), bin_edges,@view(means_matr[t,k,:]))

        res[k] = fun(@view(indiv_prop[idcs_in]))
    end
    res
end

function mean_counts_party(A::Vector{Vector{T}}, partyfori::AbstractVector) where T
    # Il numero di righe nell'array A
    n_rows = length(A)
    # Il numero di bin sarà il numero di bordi - 1
    #n_bins = length(bin_edges) - 1
    unique_part = unique(partyfori)
    nkeys = length(unique_part)

    # Pre-alloca la matrice dei risultati per le performance
    # Risultato è una matrice n_rows x n_bins di Int
    means_matr = zeros(n_rows, nkeys)
    
    # 1. Definisce l'oggetto istogramma con i bordi desiderati
    # closed=:left significa che i bin sono intervalli chiusi a sinistra: [a, b)
    #h = Histogram(bin_edges, closed=:left)
    whichids = [partyfori .== g for g in unique_part]
    
    # 2. Itera su ogni riga dell'array 2D
    for (t, row) in enumerate(A)
        # .weights contiene i conteggi per ogni bin.
        for k =1:nkeys
            #tid = tileData.tiles_idcs[k]
            #idcs_in = tileData.idcs_in_tile[tid]
            ma = whichids[k]

            means_matr[t,k] = _mean(@view(row[ma]))
        end
    end
    
    return means_matr, unique_part
end