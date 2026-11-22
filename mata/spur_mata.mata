version 16.0
mata

void spur_mata_test(string scalar flavor, string scalar yvar, string scalar xvars,
    string scalar touse, string scalar latlong, real scalar q, real scalar nrep,
    string scalar pmat, string scalar statscal, string scalar hascal, string scalar cvmat)
{
    real matrix Y, X, distmat, emat, sigma
    struct spur_diagnostic_design scalar design
    struct spur_diagnostic_result scalar result

    emat = rnormal(2*q, nrep, 0, 1)
    Y = st_data(., yvar, touse)
    distmat = getdistmat_normalized(get_s_matrix(touse, latlong))
    if (flavor == "i0resid" | flavor == "i1resid") {
        X = J(rows(Y), 1, 1)
        if (xvars != "") X = X, st_data(., xvars, touse)
        design = spur_prepare_residual_diagnostic(Y, X, distmat, q)
    }
    else {
        sigma = get_sigma_lbm_dm(distmat)
        design = spur_prepare_diagnostic(Y, sigma, get_R(sigma, q), J(0, 0, .), 1)
    }
    if (flavor == "i0" | flavor == "i0resid") result = spatial_i0_test(design, distmat, emat)
    else result = spatial_i1_test(design, distmat, emat)

    st_numscalar(pmat, result.pvalue[1])
    st_numscalar(statscal, result.LR[1])
    st_numscalar(hascal, result.ha_parm)
    st_matrix(cvmat, result.cv)
}

void spur_mata_pipeline(string scalar invars, string scalar outvars, string scalar touse,
    string scalar coordvars, real scalar latlong, real scalar q, real scalar nrep,
    string scalar diagmat, string scalar cvmat)
{
    external real scalar latlongflag
    real matrix values, Y, X, distmat, sigma, V, emat, transformed
    real colvector d
    struct spur_diagnostic_design scalar design
    struct spur_diagnostic_result scalar i0, i1, i0r, i1r

    latlongflag = latlong
    values = st_data(., invars, touse)
    Y = values[., 1]
    X = J(rows(values), 1, 1), values[., 2..cols(values)]
    emat = rnormal(2*q, nrep, 0, 1)
    distmat = getdistmat_normalized(st_data(., coordvars, touse))
    sigma = get_sigma_lbm_dm(distmat)
    spur_sorted_eigensystem(sigma, V, d)
    design = spur_prepare_diagnostic(Y, sigma, spur_select_basis(V, d, q), J(0, 0, .), 1)
    sigma = J(0, 0, .)
    transformed = spur_apply_lbm_gls(V, d, values)
    V = J(0, 0, .)
    i0 = spatial_i0_test(design, distmat, emat)
    i1 = spatial_i1_test(design, distmat, emat)

    design = spur_prepare_residual_diagnostic(Y, X, distmat, q)
    i0r = spatial_i0_test(design, distmat, emat)
    i1r = spatial_i1_test(design, distmat, emat)

    st_matrix(diagmat, (i0.LR, i0.pvalue, i0.ha_parm) \ (i1.LR, i1.pvalue, i1.ha_parm) \
        (i0r.LR, i0r.pvalue, i0r.ha_parm) \ (i1r.LR, i1r.pvalue, i1r.ha_parm))
    st_matrix(cvmat, i0.cv' \ i1.cv' \ i0r.cv' \ i1r.cv')
    st_store(., tokens(outvars), touse, transformed)
}
end
