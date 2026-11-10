*! version 1.0.2 10nov2026
version 14
mata

real scalar getcbar(real scalar avc0, real matrix distmat)
{
	real scalar denom
	real scalar dbar, floor_corr, lo, hi, c, f, slope, candidate, iter
	real matrix E
	real vector d
	d = lvech(distmat)
	denom = length(d)
	dbar = mean(d)
	floor_corr = mean(d :== 0)
	if (dbar <= 0 | dbar >= . | avc0 <= floor_corr | avc0 >= 1) {
		errprintf("Average correlation is unattainable for these locations.\n")
		_error(3498)
	}
	lo = 0
	hi = -log(avc0) / dbar
	for (iter = 1; iter <= 1024; iter++) {
		E = exp(-hi * d)
		if (mean(E) <= avc0) break
		lo = hi
		hi = 2 * hi
		if (hi >= .) _error(3498)
	}
	c = hi
	for (iter = 1; iter <= 100; iter++) {
		E = exp(-c * d)
		f = mean(E) - avc0
		if (abs(f) <= 1e-10 * avc0) return(c)
		if (f > 0) lo = c
		else hi = c
		slope = -mean(d :* E)
		candidate = c - f / slope
		if (candidate <= lo | candidate >= hi | candidate >= .) candidate = (lo + hi) / 2
		c = candidate
	}
	errprintf("Spatial scale solver did not converge.\n")
	_error(3498)
	return(.)
}

end
