*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real matrix make_transform(real matrix s, string scalar transformation, | real scalar radius, real vector cluster)
{
	real matrix H
	
	if (transformation == "nn") {
		H = nn_matrix(s)
	} else if (transformation == "iso") {
		H = iso_matrix(s, radius)
	} else if (transformation == "cluster") {
		H = cluster_matrix(cluster)
	} else {
		_error("Invalid transformation.")
	}
	
	return(H)
}

end
