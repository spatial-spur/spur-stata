*! version 1.1.0  17nov2026 spatial-spur
version 16.0
mata

void spur_sorted_eigensystem(real matrix sigma, real matrix V, real colvector d)
{
    real colvector order_by_eval
    symeigensystem(sigma, V, d)
    d = d'
    order_by_eval = order(d, -1)
    d = d[order_by_eval]
    V = V[., order_by_eval]
}

real matrix spur_select_basis(real matrix V, real colvector d, real scalar qmax)
{
    real matrix R
    real colvector ds, colj
    real scalar j, i, best_abs, tolerance

    if (qmax < 2 | qmax >= rows(V)) {
        errprintf("q must satisfy 2 <= q < n\n")
        exit(198)
    }

	ds = d[1..qmax]
	if (ds[1] <= 0 | ds[1] >= . | ds[qmax] >= . | ds[qmax] <= 1e-12 * ds[1]) {
		errprintf("SPUR covariance basis has insufficient positive rank for requested q\n")
		exit(198)
	}
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

real matrix get_R(real matrix sigma, real scalar q)
{
    real matrix V
    real colvector d
    spur_sorted_eigensystem(sigma, V, d)
    return(spur_select_basis(V, d, q))
}
end
