*! version 1.0.0  21jan2025 pdavidboll
version 16.0
mata

real scalar get_ha_parm_I0(real matrix om_i0, real matrix om_bm, real matrix e, | real scalar missed_power) {
    real matrix mixture, limit_shape
    real scalar best_param, best_power
    // Variable declarations
    real scalar g, pow, gtry, g1, g2
    real scalar ii

    // Initial search for lower bound where power < 0.5
    pow = 1
    gtry = 1
    ii = 0
    while (pow > 0.5) {
        if (++ii > 200) {
            errprintf("SPUR power calibration could not bracket the target\n")
            exit(198)
        }
        g = gtry;
        // Verify that power is less than 50%
        pow = getpow_qf(om_i0, om_i0 + g * om_bm, e)
        gtry = g / 2
    }
    missed_power = .
    if (pow == 0.5) return(g)
    g1 = g

    // Initial search for upper bound where power > 0.5
    best_param = g
    best_power = pow
    limit_shape = om_bm / trace(om_bm)
    pow = 0
    gtry = (ii > 1 ? 2*g1 : 30)
    ii = 0
    while (pow < 0.5) {
        if (++ii > 200) {
            missed_power = best_power
            return(best_param)
        }
        g = gtry
        // Verify that power is greater than 50%
        mixture = om_i0 + g * om_bm
        pow = getpow_qf(om_i0, mixture, e)
        if (pow > best_power) {
            best_power = pow
            best_param = g
        }
        if (pow < 0.5 & max(abs(mixture / trace(mixture) - limit_shape)) <= 8*epsilon(1)) {
            missed_power = best_power
            return(best_param)
        }
        gtry = g * 2
    }
    g2 = g

    // Binary search for precise value where power is approximately 0.5
    ii = 1
    while (abs(pow - 0.5) > 0.01) {
        g = (g1 + g2) / 2
        pow = getpow_qf(om_i0, om_i0 + g * om_bm, e)
        if (pow > 0.5) {
            g2 = g
        } else if (pow < 0.5) {
            g1 = g
        }
        ii = ii + 1
        if (ii > 20) {
            break
        }
    }

    return(g)
}

end