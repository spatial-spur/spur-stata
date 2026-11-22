*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real scalar getcbar(real scalar rhobar, real matrix distmat)
{
    real colvector distances
    real scalar lo, hi, mid, legacy

    distances = lvech(distmat)
    if (rhobar <= mean(distances :== 0) | rhobar >= 1) {
        errprintf("average correlation is unattainable for these locations\n")
        exit(198)
    }

    lo = hi = 10
    while (mean(exp(-lo * distances)) <= rhobar) {
        hi = lo
        lo = lo / 2
        if (lo == 0) {
            errprintf("SPUR decay calibration lower bracket underflowed\n")
            exit(198)
        }
    }
    while (mean(exp(-hi * distances)) >= rhobar) {
        lo = hi
        hi = 2 * hi
        if (hi >= .) {
            errprintf("SPUR decay calibration upper bracket exceeds the numeric range\n")
            exit(198)
        }
    }

    // Retain the legacy tolerance only when its search found a bracket.
    legacy = hi <= 10240
    while (hi - lo > (legacy ? .001 : 1e-12 * hi)) {
        mid = legacy ? sqrt(lo * hi) : lo + (hi - lo) / 2
        if (mid == lo | mid == hi) break
        if (mean(exp(-mid * distances)) < rhobar) hi = mid
        else lo = mid
    }
    return(legacy ? sqrt(lo * hi) : lo + (hi - lo) / 2)
}

end
