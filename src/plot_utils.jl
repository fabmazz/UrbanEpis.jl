function  make_matrix_heatmap(dfplot::DataFrame,i,j,value_col;reverse_i=false, reverse_j=false)

    imin = minimum(dfplot[!,i])
    jmin = minimum(dfplot[!,j])
    imax = maximum(dfplot[!,i])
    jmax = maximum(dfplot[!,j])
    ni= imax-imin+1
    nj = jmax -jmin+1
    znan::Matrix{Float64} = fill(NaN, (ni,nj))
    
    for (x,y,t) in eachrow(dfplot[!,[i,j,value_col]])
        ii = if reverse_i 
            ni -(x-imin)  
        else x-imin+1 end 

        if reverse_j
            znan[ii,end-(y-jmin)] = t
        else
            znan[ii, y-jmin +1] = t
        end
    end 

    znan
end
