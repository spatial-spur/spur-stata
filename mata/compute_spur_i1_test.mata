*! version 1.0.0  21jan2025 pdavidboll
version 16.0
mata

struct spur_diagnostic_result scalar spatial_i1_test(struct spur_diagnostic_design scalar design,
    real matrix distmat, real matrix emat)
{
    real scalar ha_parm, c01, c95, i
    real colvector LR, pvalue, P
    real rowvector lr_draws
    real matrix R, om_ho, om_ha, omi_ho, omi_ha, samples
    struct spur_diagnostic_result scalar result

    R = design.R
    om_ho = design.om_bm
    c01 = getcbar(.01, distmat)
    c95 = getcbar(.95, distmat)
    if (spur_shape_equal(om_ho, spur_diagnostic_omega(distmat, c01, R, design.M)) &
        spur_shape_equal(om_ho, spur_diagnostic_omega(distmat, getcbar(.5, distmat), R, design.M)) &
        spur_shape_equal(om_ho, spur_diagnostic_omega(distmat, c95, R, design.M))) {
        errprintf("SPUR covariance shapes are uninformative\n")
        exit(198)
    }

    ha_parm = get_ha_parm_I1(om_ho, distmat, R, emat, design.M, c01, c95, om_ha)
    omi_ho = invsym(om_ho)
    omi_ha = invsym(om_ha)
    samples = cholesky(om_ho) * emat
    lr_draws = colsum((cholesky(omi_ho)' * samples):^2) :/ colsum((cholesky(omi_ha)' * samples):^2)

    LR = J(cols(design.projected_y), 1, .)
    pvalue = J(rows(LR), 1, .)
    for (i = 1; i <= rows(LR); i++) {
        P = design.projected_y[., i]
        LR[i] = (P' * omi_ho * P) / (P' * omi_ha * P)
        pvalue[i] = mean(lr_draws' :> LR[i])
    }

    result.LR = LR
    result.pvalue = pvalue
    result.cv = spur_quantile(lr_draws', (.99 \ .95 \ .90))
    result.ha_parm = ha_parm
    return(result)
}
end
