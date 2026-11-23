*! version 1.0.0  21jan2025 pdavidboll
version 16.0
mata

struct spur_diagnostic_result scalar spatial_i0_test(struct spur_diagnostic_design scalar design,
    real matrix distmat, real matrix emat)
{
    real scalar ha_parm, i, j, q, c, floor_corr, first
    real colvector rho_grid, LR, pvalue, cv, P
    real rowvector lr_draws
    real matrix R, om_i0, om_ha, ch_omi_ho, ch_omi_ha, samples
    struct spur_diagnostic_result scalar result

    R = design.R
    q = cols(R)
    floor_corr = mean(lvech(distmat) :== 0)
    if (floor_corr > .03) {
        errprintf("I(0) null requires a co-located-pair fraction at or below 0.03. %s\n", spur_site_summary(distmat))
        exit(198)
    }
    c = (floor_corr < .001 ? getcbar(.001, distmat) : .)
    om_i0 = spur_diagnostic_omega(distmat, c, R, design.M)
    if (spur_shape_equal(om_i0, design.om_bm)) {
        errprintf("SPUR covariance shapes are uninformative\n")
        exit(198)
    }

    ha_parm = get_ha_parm_I0(om_i0, design.om_bm, emat[(q+1)..(2*q), .], result.best_power)
    om_ha = om_i0 + ha_parm * design.om_bm
    ch_omi_ho = lusolve(cholesky(om_i0), I(q))
    ch_omi_ha = lusolve(cholesky(om_ha), I(q))

    LR = J(cols(design.projected_y), 1, .)
    for (i = 1; i <= rows(LR); i++) {
        P = design.projected_y[., i]
        LR[i] = colsum((ch_omi_ho * P):^2) / colsum((ch_omi_ha * P):^2)
    }

    // Maximize each tail area and critical value over the I(0) null grid.
    first = (floor_corr < .0001 ? .0001 : floor_corr + (.03 - floor_corr) / 30)
    rho_grid = rangen(first, .03, 30)
    rho_grid = floor_corr \ select(rho_grid, rho_grid :> floor_corr)
    cv = J(3, 1, 0)
    pvalue = J(rows(LR), 1, 0)
    for (i = 1; i <= rows(rho_grid); i++) {
        c = (i == 1 ? . : getcbar(rho_grid[i], distmat))
        samples = cholesky(spur_diagnostic_omega(distmat, c, R, design.M)) * emat[1..q, .]
        lr_draws = colsum((ch_omi_ho * samples):^2) :/ colsum((ch_omi_ha * samples):^2)
        if (i == 1) cv = spur_rank_cutoffs(lr_draws')
        else cv = rowmax((cv, spur_rank_cutoffs(lr_draws')))
        for (j = 1; j <= rows(LR); j++) {
            pvalue[j] = max((pvalue[j], (1 + sum(lr_draws :>= LR[j])) / (cols(emat) + 1)))
        }
    }

    result.LR = LR
    result.pvalue = pvalue
    result.cv = cv
    result.ha_parm = ha_parm
    return(result)
}
end
