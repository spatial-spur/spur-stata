*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real matrix getdistmat(real matrix s)
// computes matrix of distances from locations
// if external latlongflag=0, distances computed as norm, otherwise latitude / longitude and haversine formula for sphere with radius 1/Pi
{	
	external real scalar latlongflag
	real matrix mat
	real scalar n,i
	real matrix d, radians
	real colvector scale, rad, cos_lat, lat_sin, lon_sin
	n=rows(s)
	mat=J(n,n,0)
	if(latlongflag==0){
		for(i=1;i<=n;i++){
			d = s :- s[i,.]
			if (hasmissing(d)) {
				errprintf("pairwise coordinate differences exceed the numeric range\n")
				exit(198)
			}
			scale = rowmax(abs(d))
			scale = scale :+ (scale :== 0)
			mat[.,i] = scale :* sqrt(rowsum((d :/ scale):^2))
		}
	}
	else {
		radians = s * (pi() / 180)
		cos_lat = cos(radians[., 1])
		for (i = 1; i <= n; i++) {
			lat_sin = sin((radians[i, 1] :- radians[., 1]) * .5)
			lon_sin = sin((radians[i, 2] :- radians[., 2]) * .5)
			rad = lat_sin :* lat_sin + ((cos_lat[i] :* cos_lat) :* lon_sin) :* lon_sin
			rad = rad :* (rad :> 0)
			rad = rad :* (rad :<= 1) :+ (rad :> 1)
			mat[., i] = (1 / pi()) * asin(sqrt(rad))
		}
	}
	if (hasmissing(mat)) {
		errprintf("pairwise distances exceed the numeric range\n")
		exit(198)
	}
	return(mat)
}	

end