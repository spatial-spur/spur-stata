version 16.0
mata

struct spur_diagnostic_design {
    real matrix R, om_bm, projected_y, M
    real rowvector info
}

struct spur_diagnostic_result {
    real colvector LR, pvalue, cv
    real scalar ha_parm, best_power
}

real matrix spur_diagnostic_omega(real matrix distmat, real scalar c, real matrix R, real matrix M)
{
    real matrix kernel, projected, omega
    projected = (rows(M) == 0 ? R :- mean(R) : M * R)
    if (missing(c)) kernel = distmat :== 0
    else if (c < .01) kernel = expm1(-c * distmat)
    else kernel = exp(-c * distmat)
    omega = projected' * (kernel * projected)
    return(0.5 * (omega + omega'))
}

struct spur_diagnostic_design scalar spur_prepare_diagnostic(real matrix Y, real matrix sigma,
    real matrix R, real matrix Qx, real colvector constant_residual)
{
    real colvector y, projected, correction
    real scalar j, anchor, input_norm, scale, allowance
    struct spur_diagnostic_design scalar design

    design.R = R
    design.om_bm = R' * sigma * R
    design.M = J(0, 0, .)
    design.projected_y = J(cols(R), cols(Y), .)
    for (j = 1; j <= cols(Y); j++) {
        y = Y[., j]
        if (min(y) == max(y)) {
            errprintf("response Y is constant; SPUR inference requires variation in Y\n")
            exit(198)
        }
        anchor = 0.5 * min(y) + 0.5 * max(y)
        y = y :- anchor
        scale = max(abs(y))
        allowance = 0
        if (rows(constant_residual) == 0) y = y :/ scale
        else {
            correction = anchor * constant_residual
            scale = max((scale, max(abs(correction))))
            y = y :/ scale + correction :/ scale
            allowance = (8 * epsilon(1) * abs(anchor) / scale) * sqrt(rows(Y))
        }
        if (cols(Qx) == 0) y = y :- mean(y)
        else {
            input_norm = norm(y)
            y = y - Qx * (Qx' * y)
            y = y - Qx * (Qx' * y)
            if (norm(y) <= 1e-12 * input_norm + allowance) {
                errprintf("residualized response Y has zero variation orthogonal to X\n")
                exit(198)
            }
        }
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
    real scalar n, k, j, pass, prev, scale, cnorm, constant_norm_sq
    real matrix Qx, M, sigma, U, Vt
    real colvector singular
    real rowvector info
    real colvector column, constant
    struct spur_diagnostic_design scalar design

    n = rows(distmat)
    k = cols(X)
    if (k == 0 | k >= n) {
        errprintf("invalid residual test design matrix dimensions\n")
        exit(198)
    }
    if (q > n - k) {
        errprintf("q=%g cannot exceed residual degrees of freedom (n - rank(X) = %g)\n", q, n-k)
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

    constant = J(n, 1, 1 / sqrt(n))
    for (pass = 1; pass <= 2; pass++) {
        constant = constant - Qx * (Qx' * constant)
    }
    constant_norm_sq = sum(constant :^ 2)
    if (constant_norm_sq > (8 * epsilon(1))^2) {
        svd(X :/ colmax(abs(X)), U, singular, Vt)
        if (sqrt(constant_norm_sq) > 64 * epsilon(1) *
            (1 + singular[1] * sqrt(sum((colsum(U)' :/ (sqrt(n) * singular)):^2)))) {
            errprintf("SPUR residual tests require a design that spans a constant\n")
            exit(198)
        }
    }

    M = I(n) - Qx * Qx'
    sigma = M * (-0.5 * distmat) * M
    constant = (constant_norm_sq <= (8 * epsilon(1))^2 ? J(0, 1, .) : constant * sqrt(n))
    design = spur_prepare_diagnostic(Y, sigma, get_R(sigma, q, info, distmat), Qx, constant)
    design.info = info
    design.M = M
    return(design)
}
end
