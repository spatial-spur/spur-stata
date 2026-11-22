*! version 1.0.0  21jan2025 pdavidboll
version 16.0
mata

struct spur_diagnostic_result scalar spatial_i0_test(struct spur_diagnostic_design scalar design,
    real matrix distmat, real matrix emat)
{
    real scalar ha_parm, i, j
    real colvector rho_grid, LR, pvalue, cv, P
    real rowvector frequencies, lr_draws
    real matrix R, om_i0, om_ha, ch_omi_ho, ch_omi_ha, samples
    struct spur_diagnostic_result scalar result

    R = design.R
    om_i0 = spur_diagnostic_omega(distmat, getcbar(.001, distmat), R, design.M)
    if (spur_shape_equal(om_i0, design.om_bm)) {
        errprintf("SPUR covariance shapes are uninformative\n")
        exit(198)
    }

    ha_parm = get_ha_parm_I0(om_i0, design.om_bm, emat)
    om_ha = om_i0 + ha_parm * design.om_bm
    ch_omi_ho = cholesky(invsym(om_i0))'
    ch_omi_ha = cholesky(invsym(om_ha))'

    LR = J(cols(design.projected_y), 1, .)
    for (i = 1; i <= rows(LR); i++) {
        P = design.projected_y[., i]
        LR[i] = colsum((ch_omi_ho * P):^2) / colsum((ch_omi_ha * P):^2)
    }

    // Maximize each tail area and critical value over the I(0) null grid.
    rho_grid = rangen(.0001, .03, 30)
    frequencies = (.99, .95, .90)
    cv = J(3, 1, 0)
    pvalue = J(rows(LR), 1, 0)
    for (i = 1; i <= rows(rho_grid); i++) {
        samples = cholesky(spur_diagnostic_omega(distmat, getcbar(rho_grid[i], distmat), R, design.M)) * emat
        lr_draws = colsum((ch_omi_ho * samples):^2) :/ colsum((ch_omi_ha * samples):^2)
        cv = rowmax((cv, spur_quantile(lr_draws', frequencies')))
        for (j = 1; j <= rows(LR); j++) {
            pvalue[j] = max((pvalue[j], mean(lr_draws' :> LR[j])))
        }
    }

    result.LR = LR
    result.pvalue = pvalue
    result.cv = cv
    result.ha_parm = ha_parm
    return(result)
}
end
