function df_fields_op(df1::DataFrame, df2::DataFrame, col_index, operation::Function)
    n1 = Set(names(df1))
    n2 = Set(names(df2))
    nkeep=intersect(n1,n2)
    if !(col_index in nkeep)
        throw(ArgumentError("Column '$col_index' not found in at least one dataframe"))
    end
    dfnew = select(df1, collect(nkeep))
    delete!(nkeep, col_index)
    for col in nkeep
        dfnew[!,col]= operation.(dfnew[!,col], df2[!,col])
    end

    dfnew
end

subtract_df_fields(df1,df2, col_index) = df_fields_op(df1,df2,col_index, -)