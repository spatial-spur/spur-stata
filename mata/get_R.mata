*! version 1.1.0  17nov2026 spatial-spur
version 14
mata

real matrix get_R(real matrix sig, real scalar qmax)
{

	real matrix V, R
	real colvector d, ds, colj
	real colvector ii
	real scalar j, i, best_abs, sign_cutoff

	if (qmax < 2 | qmax >= rows(sig)) {
		errprintf("q must satisfy 2 <= q < n\n")
		exit(198)
	}
	
	symeigensystem(sig, V, d)
	d = d'
	
	ii = order(d,-1)
	d = d[ii]
	V = V[.,ii]
	
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
		sign_cutoff = best_abs * (1 - 64 * epsilon(1))
		for (i = 1; abs(colj[i]) < sign_cutoff; i++) {}
		if (colj[i] < 0) {
			R[., j] = -colj
		}
	}

	return(R)
}

end