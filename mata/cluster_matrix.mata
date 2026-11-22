*! version 1.0.0  21jan2025 pdavidboll
version 16.0
mata 

real matrix cluster_matrix(real colvector cluster)
{
// Compute the cluster matrix 
// Input:
//   cluster: vector of cluster ids
	
	real scalar n
	real matrix clust_mat
	
	n = rows(cluster)
	
	stata("disp as text"+char(34)+"Number of observations: "+strofreal(n, "%9.0g")+", number of clusters "+strofreal(rows(uniqrows(cluster)), "%9.0g"))
		
	clust_mat = J(1,n,cluster) :== J(n,1,cluster')
	
	clust_mat = clust_mat :/ rowsum(clust_mat)
				
	clust_mat = I(n) - clust_mat
	
	return(clust_mat)
	
}

real matrix spur_demean_clusters(real matrix values, real colvector cluster)
{
    real colvector permutation, sample_rows
    real matrix panels, block
    real rowvector lo, hi, anchor, scale
    real scalar g

    permutation = order(cluster, 1)
    panels = panelsetup(cluster[permutation], 1)
    for (g = 1; g <= rows(panels); g++) {
        sample_rows = permutation[panels[g, 1]..panels[g, 2]]
        block = values[sample_rows, .]
        lo = colmin(block)
        hi = colmax(block)
        anchor = 0.5 * lo + 0.5 * hi
        scale = colmax((hi - anchor) \ (anchor - lo))
        scale = scale + (scale :== 0)
        block = (block :- anchor) :/ scale
        values[sample_rows, .] = (block :- mean(block)) :* scale
    }
    if (hasmissing(values)) {
        errprintf("Cluster demeaning produced a non-finite entry.\n")
        exit(198)
    }
    return(values)
}

end