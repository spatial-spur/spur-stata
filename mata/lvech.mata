*! version 1.0.1  26mar2026 pdavidboll
version 14
mata

real vector lvech(real matrix S) {
    real scalar n, j, pos, len
    real vector v
    
    n = rows(S)
    v = J((n * (n - 1)) / 2, 1, .)
    pos = 1

    for (j = 1; j < n; j++) {
        len = n - j
        v[pos..(pos + len - 1)] = S[(j+1)..n, j]
        pos = pos + len
    }
    
    return(v)
}

end
