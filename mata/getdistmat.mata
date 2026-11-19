*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real matrix getdistmat(real matrix s)
// computes matrix of distances from locations
// if external latlongflag=0, distances computed as norm, otherwise latitude / longitude and haversine formula for sphere with radius 1/Pi
{	
	external real scalar latlongflag
	real matrix mat
	real scalar n,i,c
	real matrix d
	real colvector scale, rad
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
	else
		{
		c=3.14159265359/180
		for(i=1;i<=n;i++){
			d=(.5*c)*(s:-s[i,.])
			rad = sin(d[.,1]):^2 + cos(c*s[i,1])*(cos(c*s[.,1]):*(sin(d[.,2]):^2))
			rad = rad :* (rad :> 0)
			rad = rad :* (rad :<= 1) :+ (rad :> 1)
			mat[.,i] = asin(sqrt(rad))/3.14159265359
		}
	}	
	if (hasmissing(mat)) {
		errprintf("pairwise distances exceed the numeric range\n")
		exit(198)
	}
	return(mat)
}	

end