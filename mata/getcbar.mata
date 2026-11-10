*! version 1.0.1  26mar2026 pdavidboll
version 14
mata

real scalar getcbar(real scalar rhobar, real matrix distmat) {
	real scalar c0, c1, v, i1, jj, cm, cbar, n, denom
	
	c0 = 10
	c1 = 10
	n = rows(distmat)
	denom = n * (n - 1)

    // Adjust c0 and c1 bounds so c0 yields v > rhobar and c1 yields v < rhobar
    i1 = 0
	jj = 0

    while (i1 == 0) {
        v = (sum(exp(-c0 * distmat)) - n) / denom
        i1 = (v > rhobar)
      
        if (i1 == 0) {
            c1 = c0
            c0 = c0 / 2
            jj = jj + 1
        }
        if (jj > 500) {
            _error("rhobar too large")
        }
    }

    // Verify that c1 yields a value less than rhobar
    i1 = 0
    jj = 0
    while (i1 == 0) {
        v = (sum(exp(-c1 * distmat)) - n) / denom
        i1 = (v < rhobar)
      
        if (i1 == 0) {
            c0 = c1
            c1 = 2 * c1
            jj = jj + 1
        }
        if (c1 > 10000) {
            i1 = 1
        }
        if (jj > 500) {
            _error("rhobar too small")
        }
    }

    // Bisection determination of cbar
    while ((c1 - c0) > 0.001) {
        cm = sqrt(c0 * c1)
        v = (sum(exp(-cm * distmat)) - n) / denom
      
        if (v < rhobar) {
            c1 = cm
        } else {
            c0 = cm
        }
    }

    cbar = sqrt(c0 * c1)
    return(cbar)
}

end
