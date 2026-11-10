*! version 1.0.1  26mar2026 pdavidboll
version 14
mata

real matrix get_R(real matrix sig, real scalar qmax)
{

	real matrix V, R
	real colvector d, ds
	real colvector ii
	
	symeigensystem(sig, V, d)
	d = d'
	
	ii = order(d,-1)
	d = d[ii]
	V = V[.,ii]
	
	qmax = min((rows(sig) - 1, qmax))
	if (qmax < 1) qmax = 1
	
	ds = d[1..qmax]
	R = V[., 1..qmax]
	
	return(R)
}

end
