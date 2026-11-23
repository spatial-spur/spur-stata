*! version 1.0.0  21jan2025 pdavidboll
version 16.0
mata

real matrix spur_apply_lbm_gls(real matrix V, real colvector d, real matrix Y)
{
    real colvector keep
    real matrix basis, transformed

    keep = d :> 1e-10
    basis = select(V, keep')
    transformed = basis * ((basis' * spur_demean_clusters(Y, J(rows(Y), 1, 1))) :/ sqrt(select(d, keep)))
    return(transformed :- mean(transformed))
}

real matrix spur_lbm_gls(real matrix s, real matrix Y)
{
    real matrix V
    real colvector d
    spur_sorted_eigensystem(get_sigma_lbm_dm(getdistmat_normalized(s)), V, d)
    return(spur_apply_lbm_gls(V, d, Y))
}
end
