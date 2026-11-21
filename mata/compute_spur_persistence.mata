*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

// Define the struct for Spatial_I1_Test
struct struct_spatial_persistence {
    real scalar hl_ci_lower
	real scalar hl_ci_upper
}

real colvector spur_halflife_accept(real colvector p, real scalar nrep, real scalar level) {
    // Compare the complement count with level to keep equality excluded.
    return(((nrep :- round(p :* nrep)) :/ nrep) :< level)
}

struct struct_spatial_persistence scalar spatial_persistence(real matrix Z, real matrix distmat, real matrix emat, real scalar level) {
	
	real scalar n_hl, n_hl_ha, n_hl_grid_ho, ci_l, ci_u
	real matrix hl_grid_ho, hl_grid_ha, c_grid_ho, c_grid_ha, pv_mat, ii, hl, Z_dm

	if (rows(emat) < 2) {
		errprintf("half-life inference requires q >= 2 weighted directions\n")
		exit(198)
	}
	if (max(Z) - min(Z) == 0) {
		errprintf("response y is constant; half-life inference requires variation in y\n")
		exit(198)
	}
	Z_dm = Z :- mean(Z)
	Z_dm = Z_dm :/ max(abs(Z_dm))
	
	n_hl = 100
	hl_grid_ho = rangen(0.001,1,n_hl)
	hl_grid_ho = hl_grid_ho \ rangen(1.01,3,30)
	hl_grid_ho = hl_grid_ho \ 100
	n_hl_grid_ho = length(hl_grid_ho)
	n_hl_ha = 50
	hl_grid_ha = rangen(0.001,1.0,n_hl_ha)
	
	c_grid_ho = -log(0.5) :/ hl_grid_ho
	c_grid_ha = -log(0.5) :/ hl_grid_ha
	
	pv_mat=c_ci(Z_dm,distmat,emat,c_grid_ho,c_grid_ha)
    ii = spur_halflife_accept(pv_mat, cols(emat), level)
    hl = select(hl_grid_ho,ii)
    ci_l = min(hl)
    ci_u = max(hl)
	
	
	// Struct output
     struct struct_spatial_persistence scalar SP
     SP.hl_ci_lower = ci_l
     SP.hl_ci_upper = ci_u

     return(SP)
	
	
	
}

real matrix c_ci(real matrix Y, real matrix distmat, real matrix emat, 
                 real vector c_grid_ho, real vector c_grid_ha) {

    // Variable declarations
    real scalar q, n_c_ho, n_c_ha, i, j
    real scalar rho_bm, c_bm, c, log_det_ho, log_ho_X, log_lr_X, informative, m_x, q_x
    real matrix sigdm_bm, R, om, sigdm, omi, om_ref
    real matrix ch_null, ch_omi, ch_om_mat, ch_omi_mat, log_det_mat
    real matrix ch_omi_mat_ha, log_det_mat_ha
    real matrix pv_mat, X
    real matrix e, xc, Xc, log_ho, log_ha_mat, log_ha_mat_X, log_lr, m_row, q_e

    // Initialization
    q = rows(emat)
    n_c_ho = length(c_grid_ho)
    n_c_ha = length(c_grid_ha)
    if (q < 2) {
        errprintf("half-life inference requires q >= 2 weighted directions\n")
        exit(198)
    }

    // BM covariance matrix
    rho_bm = 0.999
    c_bm = getcbar(rho_bm, distmat)
    sigdm_bm = get_sigma_dm(distmat, c_bm)

    // Construct R and eigenvectors
    R = get_R(sigdm_bm, q)
    X = R' * (Y :- mean(Y))
    if (norm(X) == 0 | norm(X) <= 1e-12 * norm(Y :- mean(Y))) {
        errprintf("projected SPUR vector has zero norm\n")
        exit(198)
    }

    om_ref = R' * get_sigma_dm(distmat, c_grid_ho[1]) * R
    informative = 0
    for (i = 2; i <= n_c_ho; i++) {
        if (!spur_shape_equal(om_ref, R' * get_sigma_dm(distmat, c_grid_ho[i]) * R)) informative = 1
    }
    for (i = 1; i <= n_c_ha; i++) {
        if (!spur_shape_equal(om_ref, R' * get_sigma_dm(distmat, c_grid_ha[i]) * R)) informative = 1
    }
    if (!informative) {
        errprintf("half-life covariance shapes are uninformative\n")
        exit(198)
    }

    // Compute Omega matrices for c_grid_ho
    ch_om_mat = asarray_create("real", 1)
    ch_omi_mat = asarray_create("real", 1)
    log_det_mat = J(n_c_ho, 1, .)

    for (i = 1; i <= n_c_ho; i++) {
        c = c_grid_ho[i]
        sigdm = get_sigma_dm(distmat, c)
        om = R' * sigdm * R
        om = 0.5 * (om + om')
        asarray(ch_om_mat, i, cholesky(om)')
        omi = invsym(om)
        ch_omi = cholesky(omi)'
        asarray(ch_omi_mat, i, ch_omi)
        log_det_mat[i, 1] = sum(ln(diagonal(ch_omi)))
    }

    // Compute Omega matrices for c_grid_ha
    ch_omi_mat_ha = asarray_create("real", 1)
    log_det_mat_ha = J(n_c_ha, 1, .)
    log_ha_mat_X = J(n_c_ha, 1, .)

    for (i = 1; i <= n_c_ha; i++) {
        c = c_grid_ha[i]
        sigdm = get_sigma_dm(distmat, c)
        om = R' * sigdm * R
        omi = invsym(0.5 * (om + om'))
        ch_omi = cholesky(omi)'
        asarray(ch_omi_mat_ha, i, ch_omi)
        log_det_mat_ha[i, 1] = sum(ln(diagonal(ch_omi)))
        Xc = ch_omi * X
        q_x = colsum(Xc:^2)
        if (q_x <= 0 | missing(q_x)) {
            errprintf("projected SPUR vector has zero quadratic form\n")
            exit(198)
        }
        log_ha_mat_X[i, 1] = log_det_mat_ha[i, 1] - 0.5 * q * ln(q_x)
    }

    // Compare simulated likelihood ratios with the observed ratio in log-space.
    pv_mat = J(n_c_ho, 1, .)

    for (i = 1; i <= n_c_ho; i++) {
        ch_null = asarray(ch_om_mat, i)
        e = ch_null' * emat
        ch_omi = asarray(ch_omi_mat, i)
        log_det_ho = log_det_mat[i, 1]
        xc = ch_omi * e
        q_e = colsum(xc:^2)'
        if (any(q_e :<= 0) | missing(q_e)) {
            errprintf("simulated SPUR vector has zero quadratic form\n")
            exit(198)
        }
        log_ho = log_det_ho :- 0.5 * q :* ln(q_e)
        log_ha_mat = J(cols(emat), n_c_ha, .)

        // Compute for data
        Xc = ch_omi * X
        q_x = colsum(Xc:^2)
        if (q_x <= 0 | missing(q_x)) {
            errprintf("projected SPUR vector has zero quadratic form\n")
            exit(198)
        }
        log_ho_X = log_det_ho - 0.5 * q * ln(q_x)

        for (j = 1; j <= n_c_ha; j++) {
            ch_omi = asarray(ch_omi_mat_ha, j)
            xc = ch_omi * e
            q_e = colsum(xc:^2)'
            if (any(q_e :<= 0) | missing(q_e)) {
                errprintf("simulated SPUR vector has zero quadratic form\n")
                exit(198)
            }
            log_ha_mat[., j] = log_det_mat_ha[j, 1] :- 0.5 * q :* ln(q_e)
        }

        m_row = rowmax(log_ha_mat)
        log_lr = m_row :+ ln(mean(exp(log_ha_mat :- m_row)')') :- log_ho
        m_x = max(log_ha_mat_X)
        log_lr_X = m_x + ln(mean(exp(log_ha_mat_X :- m_x))) - log_ho_X
        pv_mat[i, 1] = mean(log_lr :> log_lr_X)
    }

    return(pv_mat)
}
end