*! version 1.1.0  17nov2026 spatial-spur
version 16.0
mata

string scalar spur_site_summary(real matrix D)
{
    real scalar i, sites
    sites = 1
    for (i=2; i<=rows(D); i++) sites = sites + !any(D[i, 1..(i-1)] :== 0)
    return(sprintf("%g distinct locations, co-located-pair fraction %8.6f.", sites, mean(lvech(D) :== 0)))
}

void spur_sorted_eigensystem(real matrix sigma, real matrix V, real colvector d)
{
    real colvector order_by_eval
    symeigensystem(sigma, V, d)
    d = d'
    order_by_eval = order(d, -1)
    d = d[order_by_eval]
    V = V[., order_by_eval]
}

real matrix spur_select_basis(real matrix V, real colvector d, real scalar qmax, | real rowvector info, real matrix distances)
{
    real matrix R
    real colvector ds, colj
    real scalar j, i, best_abs, tolerance

    if (qmax != trunc(qmax) | qmax < 2 | qmax >= rows(V)) {
        errprintf("q=%g must satisfy 2 <= q < n=%g\n", qmax, rows(V))
        exit(198)
    }

    ds = d[1..qmax]
    tolerance = max((1.6e-11 * d[1], 8 * epsilon(1) * rows(V)))
    if (ds[1] <= 0 | ds[1] >= . | ds[qmax] >= . | ds[qmax] <= 1e-12 * ds[1]) {
        errprintf("SPUR covariance basis has insufficient positive rank for requested q=%g. Positive rank is %g.\n", qmax, sum(d :> 1e-12 * ds[1]))
        if (args() >= 5) errprintf("%s\n", spur_site_summary(distances))
        exit(198)
    }
    if (d[qmax] - d[qmax + 1] <= tolerance) {
        errprintf("Requested q splits a tied eigenspace (q=%g). Choose q at an eigenvalue gap.\n", qmax)
        if (args() >= 5) errprintf("%s\n", spur_site_summary(distances))
        exit(198)
    }
    info = (any(ds[1..(qmax-1)] - ds[2..qmax] :<= tolerance), d[qmax+1] <= 1e-12*ds[1])
    R = V[., 1..qmax]

	// Earliest nearly-maximal entry sets the sign in both backends.
	for (j = 1; j <= cols(R); j++) {
		colj = R[., j]
		best_abs = max(abs(colj))
		tolerance = 1e-12 * best_abs
		for (i = 1; best_abs - abs(colj[i]) > tolerance; i++) {}
		if (colj[i] < 0) {
			R[., j] = -colj
		}
	}

	return(R)
}

real matrix get_R(real matrix sigma, real scalar q, | real rowvector info, real matrix distances)
{
    real matrix V
    real colvector d
    spur_sorted_eigensystem(sigma, V, d)
    if (args() == 4) return(spur_select_basis(V, d, q, info, distances))
    return(spur_select_basis(V, d, q, info))
}
end
