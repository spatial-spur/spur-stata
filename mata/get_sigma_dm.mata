*! version 1.0.0  21jan2025 pdavidboll
version 16.0
mata

real matrix get_sigma_dm(real matrix distmat, real scalar c)
{
    return(demeanmat(expm1(-c * distmat)))
}

end
