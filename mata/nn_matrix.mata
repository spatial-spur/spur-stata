*! version 1.0.0  21jan2025 pdavidboll
version 16.0
mata

real matrix nn_matrix(real matrix s)
{
    real matrix distances, neighbors
    distances = editvalue(getdistmat(s), 0, .)
    if (missing(max(distances))) {
        errprintf("Nearest-neighbor transform requires at least two distinct locations\n")
        exit(198)
    }

    // Missing entries exclude zero distances from the row minima.
    neighbors = distances :== rowmin(distances)
    return(I(rows(distances)) - neighbors :/ rowsum(neighbors))
}

end
