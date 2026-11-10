*! version 1.0.1  26mar2026 pdavidboll
version 14
mata

// Define the struct for Spatial_I1_Test
struct struct_spatial_i1_test_residual {
    real vector LR
    real vector pvalue
    real vector cv_vec
    real scalar ha_parm
}

struct struct_spatial_i1_test_residual scalar spatial_i1_test_residual(real matrix Y, real matrix X_in, real matrix distmat, real matrix emat) {
    // Variable declarations
    real scalar q, n, n_y, ha_parm, i
    real matrix sigma, sigma_dm, sigdm_bm, R, RM, lam, om_ho, sigdm_ha, om_ha
    real matrix ch_om_ho, omi_ho, omi_ha, ch_omi_ho, ch_omi_ha
    real matrix y_ho, y_ho_ho, y_ho_ha
    real vector q_ho_ho, q_ho_ha, lr_ho, sz_vec, cv_vec, LR, pvalue, X, P
	real matrix M
	real scalar rho_bm, c_bm

     // Initialize scalars
     q = rows(emat)
     n = rows(distmat)
     n_y = cols(Y)

	 // Form M
	 M = I(n)-X_in*(invsym(X_in'*X_in))*X_in'

	 // BM covariance matrix (approximation for demeaned value)
	 rho_bm = 0.999
	 c_bm = getcbar(rho_bm,distmat)
	 sigdm_bm = get_sigma_residual(distmat,c_bm,M)

     // Construct eigenvectors R and low-frequency weights
     R = get_R(sigdm_bm, q)
     RM = M * R

     // Calculate matrices used in analysis
     om_ho = cross(RM, exp(-c_bm * distmat) * RM)

     // Determine Ha parameter that yields ~50% power
     ha_parm = get_ha_parm_I1_residual(om_ho, distmat, R, emat, M)
     om_ha = cross(RM, exp(-ha_parm * distmat) * RM)

     // Draws of LR under H0
     ch_om_ho = cholesky(om_ho)'
     omi_ho = invsym(om_ho)
     omi_ha = invsym(om_ha)
     ch_omi_ho = cholesky(omi_ho)'
     ch_omi_ha = cholesky(omi_ha)'
     y_ho = ch_om_ho' * emat
     y_ho_ho = ch_omi_ho * y_ho
     y_ho_ha = ch_omi_ha * y_ho
     q_ho_ho = colsum(y_ho_ho:^2)
     q_ho_ha = colsum(y_ho_ha:^2)
     lr_ho = q_ho_ho :/ q_ho_ha

     // Critical values
     sz_vec = (0.01, 0.05, 0.10)
     cv_vec = mm_quantile(lr_ho',1,(1 :- sz_vec'))

     // Realized value of LR
     LR = J(n_y, 1, .)
     pvalue = J(n_y, 1, .)
     for (i = 1; i <= n_y; i++) {
         X = Y[,i] :- mean(Y[,i])
         P = R' * X
         y_ho_ho = ch_omi_ho * P
         y_ho_ha = ch_omi_ha * P
         q_ho_ho = colsum(y_ho_ho:^2)'
         q_ho_ha = colsum(y_ho_ha:^2)'
         LR[i] = q_ho_ho / q_ho_ha
         pvalue[i] = mean(lr_ho' :> LR[i])
     }

     // Struct output
     struct struct_spatial_i1_test_residual scalar SP
     SP.LR = LR
     SP.pvalue = pvalue
     SP.cv_vec = cv_vec
     SP.ha_parm = ha_parm

     return(SP)
}

end
