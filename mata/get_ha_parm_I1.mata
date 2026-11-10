*! version 1.0.1  26mar2026 pdavidboll
version 14
mata

real scalar get_ha_parm_I1(real matrix om_ho, real matrix distmat, real matrix R, real matrix e) {
    // Variable declarations
    real scalar pow, ctry, c, pow50, c1, c2, ii
    real matrix om_c, om_c1, om_c2

    // Initial values
    pow50 = 0.5
    pow = 1
    ctry = getcbar(0.95, distmat)  // Initial guess for c

    // Step 1: Decrease `c` until `pow` is less than 0.5
    ii = 0
    while (pow > pow50) {
        ii = ii + 1
        if (ii > 40) {
            errprintf("Cannot bracket 50% power for the I(1) alternative; increase q or nrep.\n")
            _error(3498)
        }
        c = ctry
        om_c = cross(R, exp(-c * distmat) * R)
        pow = getpow_qf(om_ho, om_c, e)
        ctry = ctry / 2
    }
    c1 = c
    om_c1 = om_c

    // Step 2: Increase `c` until `pow` is greater than 0.5
    pow = 0
    ctry = getcbar(0.01, distmat)  // Another initial guess for c
    ii = 0
    while (pow < pow50) {
        ii = ii + 1
        if (ii > 40) {
            errprintf("Cannot bracket 50% power for the I(1) alternative; increase q or nrep.\n")
            _error(3498)
        }
        c = ctry
        om_c = cross(R, exp(-c * distmat) * R)
        pow = getpow_qf(om_ho, om_c, e)
        ctry = 2 * ctry
    }
    c2 = c
    om_c2 = om_c

    // Step 3: Bisection method to find `c` such that `pow` is approximately 0.5
    ii = 1
    while (abs(pow - pow50) > 0.01) {
        c = (c1 + c2) / 2
        om_c = cross(R, exp(-c * distmat) * R)
        pow = getpow_qf(om_ho, om_c, e)

        if (pow > pow50) {
            c2 = c
        } else if (pow < pow50) {
            c1 = c
        }

        ii = ii + 1
        if (ii > 20) {
            break
        }
    }

    return(c)
}

end
