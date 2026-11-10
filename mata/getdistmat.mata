*! version 1.0.1  26mar2026 pdavidboll
version 14
mata

real matrix getdistmat(real matrix s)
// computes matrix of distances from locations
// if external latlongflag=0, distances computed as norm, otherwise latitude / longitude and haversine formula for sphere with radius 1/Pi
{
	external real scalar latlongflag
	real matrix mat
	real scalar n, i, c, pi_const
	real vector lat, lon, cos_lat, d1, d2, col, h
	real matrix idx

	n = rows(s)
	mat = J(n, n, 0)

	if (latlongflag == 0) {
		for (i = 1; i < n; i++) {
			idx = (i+1)..n
			col = sqrt(rowsum((s[idx, .] :- s[i, .]):^2))
			mat[idx, i] = col
			mat[i, idx] = col'
		}
	}
	else {
		c = 3.14159265359 / 180
		pi_const = 3.14159265359
		lat = c * s[., 1]
		lon = c * s[., 2]
		cos_lat = cos(lat)
		for (i = 1; i < n; i++) {
			idx = (i+1)..n
			d1 = 0.5 * (lat[idx] :- lat[i])
			d2 = 0.5 * (lon[idx] :- lon[i])
			h = sin(d1):^2 + (cos_lat[i] * cos_lat[idx]) :* (sin(d2):^2)
            h = h :* (h :> 0)
            h = h :* (h :<= 1) :+ (h :> 1)
            col = asin(sqrt(h)) / pi_const
			mat[idx, i] = col
			mat[i, idx] = col'
		}
	}
	return(mat)
}

end
