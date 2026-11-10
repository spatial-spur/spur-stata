*! version 1.0.1  26mar2026 pdavidboll
version 14
mata 

real matrix get_sigma_lbm_dm(real matrix distmat)
{	
// Double-demeaned LBM covariance matrix:
// sigma_lbm = 0.5 * (d1*1' + 1*d1' - distmat)
// Since M1 * 1 = 0, origin terms vanish under double demeaning:
// demeanmat(sigma_lbm) == -0.5 * demeanmat(distmat)

	return(-0.5 * demeanmat(distmat))
}	
end
