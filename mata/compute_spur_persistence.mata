version 16.0
mata

struct struct_spatial_persistence {
    real scalar hl_ci_lower, hl_ci_upper
    real matrix components
}

real matrix spur_half_life_covariance(real matrix D, real matrix R, real scalar t)
{
    real scalar c
    real matrix K, omega
    if (t == 0) K = D :== 0
    else if (t == 1) K = -0.5 * D
    else {
        c = ln(2) * (1 - t) / t
        if (c < .01) K = 0.5 * expm1(-c * D) / c
        else K = exp(-ln(2) * editmissing(D / (t / (1 - t)), 1100))
    }
    omega = R' * (K * R)
    return(0.5 * (omega + omega'))
}

real scalar spur_half_life_count(real matrix omega, real matrix factors, real colvector log_weights,
    real matrix emat, real rowvector draw_norm, real colvector observed, real scalar observed_mix)
{
    real scalar q, j, threshold
    real matrix e, lower, precision, log_terms
    real rowvector peak, scores
    q = rows(emat)
    lower = cholesky(omega)
    e = lower * emat
    precision = lusolve(lower, I(q))
    threshold = observed_mix + 0.5 * q * ln(sum((precision * observed):^2))
    log_terms = J(rows(log_weights), cols(emat), .)
    for (j = 1; j <= rows(log_weights); j++) {
        precision = factors[(q * (j - 1) + 1)..(q * j), .]
        log_terms[j, .] = log_weights[j] :- 0.5 * q * ln(colsum((precision * e):^2))
    }
    peak = colmax(log_terms)
    scores = peak + ln(colsum(exp(log_terms :- peak))) + draw_norm
    return(sum(scores :>= threshold))
}

struct struct_spatial_persistence scalar spatial_persistence(real matrix Z, real matrix distmat,
    real matrix emat, real scalar level)
{
    real scalar q, nrep, i, j, k, left, right, width, h, start, mid, accepted, observed_mix, peak
    real scalar near_zero, length, tolerance
    real matrix R, bm, omega0, omega, precision, factors, bounds
    real rowvector nodes, weights, draw_norm, panels
    real colvector observed, log_weights, log_observed, grid, status
    struct spur_diagnostic_design scalar design
    struct struct_spatial_persistence scalar result

    q = rows(emat)
    nrep = cols(emat)
    bm = get_sigma_lbm_dm(distmat)
    R = get_R(bm, q)
    design = spur_prepare_diagnostic(Z, bm, R, J(0, 0, .), J(0, 1, .))
    observed = design.projected_y
    R = R :- mean(R)
    omega0 = spur_half_life_covariance(distmat, R, 0)
    if (spur_shape_equal(omega0, spur_half_life_covariance(distmat, R, 1))) {
        // If the endpoint shapes coincide, check the interior alternatives as well.
        accepted = 0
    }
    else accepted = 1

    nodes = (-0.906179845938664, -0.5384693101056831, 0, 0.5384693101056831, 0.906179845938664)
    weights = (0.2369268850561891, 0.4786286704993665, 128 / 225, 0.4786286704993665, 0.2369268850561891)
    near_zero = max((min(select(vec(distmat), vec(distmat) :> 0)) / 32, 2.2250738585072014e-308))
    right = 1 / 512
    while (right > near_zero) right = right / 2
    panels = J(1, 0, .)
    for (; right <= 1; right = 2 * right) panels = panels, right
    length = 5 * cols(panels)
    factors = J(length * q, q, .)
    log_weights = J(length, 1, .)
    log_observed = J(length, 1, .)
    left = 0
    for (i = 1; i <= cols(panels); i++) {
        right = panels[i]
        width = (right - left) / 2
        for (j = 1; j <= 5; j++) {
            k = (i - 1) * 5 + j
            h = left + width * (1 + nodes[j])
            omega = spur_half_life_covariance(distmat, R, h / (1 + h))
            if (!spur_shape_equal(omega0, omega)) accepted = 1
            precision = lusolve(cholesky(omega), I(q))
            factors[(q * (k - 1) + 1)..(q * k), .] = precision
            log_weights[k] = ln(width * weights[j]) + sum(ln(diagonal(precision)))
            log_observed[k] = log_weights[k] - 0.5 * q * ln(sum((precision * observed):^2))
        }
        left = right
    }
    if (!accepted) {
        errprintf("Half-life experiment is uninformative: covariance shapes are proportional.\n")
        exit(198)
    }
    peak = max(log_observed)
    observed_mix = peak + ln(sum(exp(log_observed :- peak)))
    draw_norm = 0.5 * q * ln(colsum(emat:^2))
    grid = rangen(1 / 128, 1, 128)
    h = 1 / 127
    while (h > near_zero) {
        h = h * 0.8408964152537145
        grid = (h / (1 + h)) \ grid
    }
    grid = 0 \ grid
    status = J(rows(grid), 1, .)
    for (i = 1; i <= rows(grid); i++) {
        omega = spur_half_life_covariance(distmat, R, grid[i])
        k = spur_half_life_count(omega, factors, log_weights, emat, draw_norm, observed, observed_mix)
        status[i] = (nrep - k) / (nrep + 1) < level
    }

    bounds = J(0, 2, .)
    start = 0
    for (i = 1; i < rows(grid); i++) {
        if (status[i] == status[i + 1]) continue
        left = grid[i]
        right = grid[i + 1]
        tolerance = min((1e-7, 1e-5 * (right - left)))
        while (right - left > tolerance) {
            mid = (left + right) / 2
            omega = spur_half_life_covariance(distmat, R, mid)
            k = spur_half_life_count(omega, factors, log_weights, emat, draw_norm, observed, observed_mix)
            accepted = (nrep - k) / (nrep + 1) < level
            if (accepted == status[i]) left = mid
            else right = mid
        }
        if (status[i]) bounds = bounds \ (start, right)
        else start = left
    }
    if (status[rows(grid)]) bounds = bounds \ (start, 1)
    result.components = bounds :/ (1 :- bounds)
    result.hl_ci_lower = rows(bounds) ? result.components[1, 1] : .
    result.hl_ci_upper = rows(bounds) ? result.components[rows(bounds), 2] : .
    return(result)
}
end
