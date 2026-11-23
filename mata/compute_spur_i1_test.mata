*! version 1.0.0  21jan2025 pdavidboll
version 16.0
mata

struct spur_diagnostic_result scalar spatial_i1_test(struct spur_diagnostic_design scalar design,
    real matrix distmat, real matrix emat)
{
    real scalar ha_parm, c01, c50, c95, i, q, floor_corr
    real colvector LR, pvalue, P
    real rowvector lr_draws
    real matrix R, om_ho, om_ha, lower, precision_ho, precision_ha, samples
    struct spur_diagnostic_result scalar result

    R = design.R
    q = cols(R)
    om_ho = design.om_bm
    floor_corr = mean(lvech(distmat) :== 0)
    c01 = getcbar(floor_corr + (1 - floor_corr) * .01, distmat)
    c50 = getcbar(floor_corr + (1 - floor_corr) * .5, distmat)
    c95 = getcbar(floor_corr + (1 - floor_corr) * .95, distmat)
    if (spur_shape_equal(om_ho, spur_diagnostic_omega(distmat, c01, R, design.M)) &
        spur_shape_equal(om_ho, spur_diagnostic_omega(distmat, c50, R, design.M)) &
        spur_shape_equal(om_ho, spur_diagnostic_omega(distmat, c95, R, design.M))) {
        errprintf("SPUR covariance shapes are uninformative\n")
        exit(198)
    }

    ha_parm = get_ha_parm_I1(om_ho, distmat, R, emat[(q+1)..(2*q), .], design.M, c01, c95, om_ha, result.best_power)
    lower = cholesky(om_ho)
    samples = lower * emat[1..q, .]
    precision_ho = lusolve(lower, I(q))
    precision_ha = lusolve(cholesky(om_ha), I(q))
    lr_draws = colsum((precision_ho * samples):^2) :/ colsum((precision_ha * samples):^2)

    LR = J(cols(design.projected_y), 1, .)
    pvalue = J(rows(LR), 1, .)
    for (i = 1; i <= rows(LR); i++) {
        P = design.projected_y[., i]
        LR[i] = sum((precision_ho * P):^2) / sum((precision_ha * P):^2)
        pvalue[i] = (1 + sum(lr_draws :>= LR[i])) / (cols(emat) + 1)
    }

    result.LR = LR
    result.pvalue = pvalue
    result.cv = spur_rank_cutoffs(lr_draws')
    result.ha_parm = ha_parm
    return(result)
}
end
