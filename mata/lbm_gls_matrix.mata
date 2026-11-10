*! version 1.0.1  26mar2026 pdavidboll
version 14
mata 

real matrix lbm_gls_matrix(real matrix s)
{
// Compute the lbm_gls matrix 
// Input:
//   s: (n,2) matrix of locations
//   latlongflag: 1 if the input is latitute and longitude, 0 if euclidean distance is to be used
	
	external real scalar latlongflag
	real matrix distmat
	real matrix sigma_lbm_dm
	real matrix V
	real colvector d
	real colvector eval
	real matrix evec
	real colvector ii
	real colvector dsi
	real matrix LBMGLS_mat
	real scalar small
		
	small = 1.0e-10
		
	distmat = getdistmat_normalized(s)
	sigma_lbm_dm = get_sigma_lbm_dm(distmat)
	
	symeigensystem(sigma_lbm_dm, V, d)
	d = d'
	
	ii = order(d,-1)
	eval = d[ii]
	evec = V[.,ii]
	
	ii = eval :> small
	eval = select(eval, ii)
	evec = select(evec, ii')
	dsi = 1 :/ sqrt(eval)
	LBMGLS_mat = (evec :* dsi') * evec'
	
	return(LBMGLS_mat)
}
end
