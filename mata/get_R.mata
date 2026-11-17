*! version 1.1.0  17nov2026 spatial-spur
version 14
mata

real matrix get_R(real matrix sig, real scalar qmax)
{

	real matrix V, R
	real colvector d, ds, colj
	real colvector ii
	real scalar j, i, best_i, best_abs, v

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

	// Canonicalize column signs (largest-magnitude entry positive, earliest index breaks ties).
	for (j = 1; j <= cols(R); j++) {
		colj = R[., j]
		best_i = 1
		best_abs = abs(colj[1])
		for (i = 2; i <= rows(colj); i++) {
			v = abs(colj[i])
			if (v > best_abs) {
				best_abs = v
				best_i = i
			}
		}
		if (colj[best_i] < 0) {
			R[., j] = -colj
		}
	}
	
	return(R)
}

end