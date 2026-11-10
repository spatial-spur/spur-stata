*! version 1.0.2  10nov2026
version 14
mata

real matrix get_R(real matrix sig, real scalar qmax)
{
    real matrix V
    real rowvector evals
    real colvector d, ii
    if (qmax < 1 | qmax != floor(qmax) | qmax >= rows(sig)) _error(3498)
    symeigensystem(sig, V, evals)
    d = evals'
    ii = order(d, -1)
    d = d[ii]
    if (d[qmax] <= max(abs(d)) * 1e-10) {
        errprintf("q exceeds the positive rank of the spatial covariance.\n")
        _error(3498)
    }
    return(V[., ii[1..qmax]])
}

end
