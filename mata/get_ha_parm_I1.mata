*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real scalar get_ha_parm_I1(real matrix om_ho, real matrix distmat, real matrix R, real matrix e,
    real matrix M, real scalar c01, real scalar c95, real matrix om_c)
{
    // Variable declarations
    real scalar pow, ctry, c, pow50, c1, c2, ii

    // Initial values
    pow50 = 0.5
    pow = 1
    ctry = c95

    // Step 1: Decrease `c` until `pow` is less than 0.5
    ii = 0
    while (pow > pow50) {
        if (++ii > 200) {
            errprintf("SPUR power calibration could not bracket the target\n")
            exit(198)
        }
        c = ctry
        om_c = spur_diagnostic_omega(distmat, c, R, M)
        pow = getpow_qf(om_ho, om_c, e)
        ctry = ctry / 2
    }
    c1 = c

    // Step 2: Increase `c` until `pow` is greater than 0.5
    pow = 0
    ctry = c01
    ii = 0
    while (pow < pow50) {
        if (++ii > 200) {
            errprintf("SPUR power calibration could not bracket the target\n")
            exit(198)
        }
        c = ctry
        om_c = spur_diagnostic_omega(distmat, c, R, M)
        pow = getpow_qf(om_ho, om_c, e)
        ctry = 2 * ctry
    }
    c2 = c

    // Step 3: Bisection method to find `c` such that `pow` is approximately 0.5
    ii = 1
    while (abs(pow - pow50) > 0.01) {
        c = (c1 + c2) / 2
        om_c = spur_diagnostic_omega(distmat, c, R, M)
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