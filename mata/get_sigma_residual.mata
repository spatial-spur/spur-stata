*! version 1.0.0  21jan2025 pdavidboll
version 16.0
mata

real matrix get_sigma_residual(real matrix distmat, real scalar c, real matrix M)
{
    real colvector constant

    // Restore the constant term for general projection matrices.
    constant = rowsum(M)
    return(M * expm1(-c * distmat) * M' + constant * constant')
}

end
