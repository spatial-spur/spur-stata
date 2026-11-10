*! version 1.0.1  26mar2026 pdavidboll
version 14
mata

real matrix demeanmat(real matrix mat)
// demeans matrix without mutating input argument
{
	real matrix res
	res = mat :- (rowsum(mat) / cols(mat))
	res = res :- (colsum(res) / rows(mat))
	return(res)
}

end
