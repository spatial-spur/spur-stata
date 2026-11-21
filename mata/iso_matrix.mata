*! version 1.0.0  21jan2025 pdavidboll
version 14
mata 

real matrix iso_matrix(real matrix s, real scalar b)
{
// Compute the iso matrix 
// Input:
//   s: (n,2) matrix of locations
//   b: radius (in metres if latlong, in units of coordinates otherwise)
//   latlongflag: 1 if the input is latitute and longitude, 0 if euclidean distance is to be used
	
	external real scalar latlongflag
	real matrix distmat
	real matrix dist_below_b
	real matrix iso_mat
	real scalar num_no_neighbors

		
	distmat = getdistmat(s) 
	
	if(latlongflag==1) {
		distmat = distmat * 2 * pi() * 6371000.009 // great circle distance
	}
	
	stata("disp as text"+char(34)+"Sanity check: Distance between observations 1 and 2 is "+strofreal(distmat[1,2], "%9.0g")+", specified radius is "+strofreal(b, "%9.0g"))
	
	dist_below_b = (distmat :> 0) :& (distmat :<= b)

	num_no_neighbors = sum(rowsum(dist_below_b) :== 0)

	if (num_no_neighbors > 0) {
		stata("disp as text"+char(34)+strofreal(num_no_neighbors, "%9.0g")+" observations have no neighbors within radius, will be set to 0.")
	}

	dist_below_b = dist_below_b :/ rowsum(dist_below_b)

	iso_mat = I(rows(distmat)) - dist_below_b

	iso_mat = editmissing(iso_mat, 0) // replace those without neighbors
	
	return(iso_mat)
	
}
end