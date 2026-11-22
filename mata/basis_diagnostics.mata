version 14
mata

struct spur_diagnostic_design {
    real matrix R, om_bm, projected_y, M
}

struct spur_diagnostic_result {
    real colvector LR, pvalue, cv
    real scalar ha_parm
}

real matrix spur_diagnostic_omega(real matrix distmat, real scalar c, real matrix R, real matrix M)
{
    real matrix sigma
    if (rows(M) == 0) sigma = get_sigma_dm(distmat, c)
    else sigma = get_sigma_residual(distmat, c, M)
    return(R' * sigma * R)
}

struct spur_diagnostic_design scalar spur_prepare_diagnostic(real matrix Y, real matrix sigma,
    real matrix R, real matrix M)
{
    real colvector y, projected
    real scalar j
    struct spur_diagnostic_design scalar design

    design.R = R
    design.om_bm = R' * sigma * R
    design.M = M
    design.projected_y = J(cols(R), cols(Y), .)
    for (j = 1; j <= cols(Y); j++) {
        y = Y[., j] :- Y[1, j]
        if (hasmissing(y)) {
            y = Y[., j] :/ max(abs(Y[., j]))
            y = y :- y[1]
        }
        if (max(abs(y)) == 0) {
            errprintf("response Y is constant; SPUR inference requires variation in Y\n")
            exit(198)
        }
        y = y :/ max(abs(y))
        y = y :- mean(y)
        projected = R' * y
        if (norm(projected) <= 1e-12 * norm(y)) {
            errprintf("projected SPUR vector has zero norm\n")
            exit(198)
        }
        design.projected_y[., j] = projected
    }
    return(design)
}

struct spur_diagnostic_design scalar spur_prepare_residual_diagnostic(real matrix Y,
    real matrix X, real matrix distmat, real scalar q)
{
    real scalar n, k, j, pass, prev, scale, cnorm
    real matrix Qx, M, sigma
    real colvector column

    n = rows(distmat)
    k = cols(X)
    if (k == 0 | k >= n) {
        errprintf("invalid residual test design matrix dimensions\n")
        exit(198)
    }
    if (q > n - k) {
        errprintf("q cannot exceed residual degrees of freedom (n - rank(X))\n")
        exit(198)
    }

    // Column-normalized MGS with reorthogonalization.
    Qx = J(n, k, 0)
    for (j = 1; j <= k; j++) {
        column = X[., j]
        scale = max(abs(column))
        if (scale <= 0 | scale >= .) {
            errprintf("design matrix is rank-deficient; residual tests require full column rank\n")
            exit(198)
        }
        column = column :/ scale
        Qx[., j] = column :/ sqrt(colsum(column :^ 2))
    }
    for (j = 1; j <= k; j++) {
        for (pass = 1; pass <= 2; pass++) {
            for (prev = 1; prev < j; prev++) {
                Qx[., j] = Qx[., j] - Qx[., prev] * (Qx[., prev]' * Qx[., j])
            }
        }
        cnorm = sqrt(colsum(Qx[., j] :^ 2))
        if (cnorm >= . | cnorm <= 1e-10) {
            errprintf("design matrix is rank-deficient; residual tests require full column rank\n")
            exit(198)
        }
        Qx[., j] = Qx[., j] :/ cnorm
    }

    M = I(n) - Qx * Qx'
    // The intercept projection removes the all-ones term from exp(-cD).
    sigma = M * expm1(-getcbar(.999, distmat) * distmat) * M
    return(spur_prepare_diagnostic(Y, sigma, get_R(sigma, q), M))
}
end
