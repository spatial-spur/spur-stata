*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real scalar get_ha_parm_I0(real matrix om_ho, real matrix om_i0, real matrix om_bm, real matrix e) {
    // Variable declarations
    real scalar g, pow, gtry, g1, g2
    real scalar ii
    real scalar si, sb, q
    real matrix L, Wi, Wb

    // Proportional covariances carry no information once variance is profiled out.
    L = cholesky(om_ho)
    Wi = lusolve(L, om_i0)
    Wi = lusolve(L, Wi')'
    Wb = lusolve(L, om_bm)
    Wb = lusolve(L, Wb')'
    q = rows(om_ho)
    si = trace(Wi) / q
    sb = trace(Wb) / q
    if (si > 0 & sb > 0 & max(abs(Wi / si - I(q))) <= 1e-10 & max(abs(Wb / sb - I(q))) <= 1e-10) {
        errprintf("Cannot bracket 50% power: the I(0) alternative has no covariance-shape information.\n")
        _error(3498)
    }

    // Initial search for lower bound where power < 0.5
    pow = 1
    gtry = 1
    ii = 0
    while (pow > 0.5) {
        ii = ii + 1
        if (ii > 64) {
            errprintf("Cannot bracket 50% power for the I(0) alternative; increase q or nrep.\n")
            _error(3498)
        }
        g = gtry;
        // Verify that power is less than 50%
        pow = getpow_qf(om_ho, om_i0 + g * om_bm, e)
        if (pow >= .) _error(3498)
        gtry = g / 2
    }
    g1 = g;

    // Initial search for upper bound where power > 0.5
    pow = 0
    gtry = 30
    ii = 0
    while (pow < 0.5) {
        ii = ii + 1
        if (ii > 64) {
            errprintf("Cannot bracket 50% power for the I(0) alternative; increase q or nrep.\n")
            _error(3498)
        }
        g = gtry
        // Verify that power is greater than 50%
        pow = getpow_qf(om_ho, om_i0 + g * om_bm, e)
        if (pow >= .) _error(3498)
        gtry = g * 2
    }
    g2 = g

    // Binary search for precise value where power is approximately 0.5
    ii = 1
    while (abs(pow - 0.5) > 0.01) {
        g = (g1 + g2) / 2
        pow = getpow_qf(om_ho, om_i0 + g * om_bm, e)
        if (pow >= .) _error(3498)
        if (pow > 0.5) {
            g2 = g
        } else if (pow < 0.5) {
            g1 = g
        }
        ii = ii + 1
        if (ii > 60 & abs(pow - 0.5) > 0.01) {
            errprintf("Cannot calibrate the I(0) alternative to 50% power; increase q or nrep.\n")
            _error(3498)
        }
    }

    return(g)
}

end
