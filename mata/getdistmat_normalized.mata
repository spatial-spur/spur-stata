*! version 1.0.0  21jan2025 pdavidboll
version 14
mata 

real matrix getdistmat_normalized(real matrix s)
{	
// Get matrix of distances from location matrix (using scpc function)
// Normalize the matrix so the largest distance is 1

	external real scalar latlongflag
	real matrix distmat
	real scalar max_d
	
	distmat = getdistmat(s)
	max_d = max(distmat)
	if (max_d <= 0 | max_d >= .) {
		errprintf("all coordinates are identical; normalization requires distinct locations\n")
		exit(198)
	}
	distmat = distmat :/ max_d
	
	return(distmat)
}	
end