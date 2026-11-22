version 16.0
mata
real colvector spur_quantile(real colvector x, real vector probs)
{
    real colvector sorted, out
    real scalar i, pos, lo, hi
    sorted = sort(x, 1)
    out = J(length(probs), 1, .)
    for (i = 1; i <= length(probs); i++) {
        pos = probs[i] * (rows(sorted) - 1) + 1
        lo = floor(pos)
        hi = ceil(pos)
        out[i] = sorted[lo] + (pos - lo) * (sorted[hi] - sorted[lo])
    }
    return(out)
}

real colvector spur_rank_cutoffs(real colvector draws)
{
    real colvector rank, cutoffs
    real scalar j
    rank = floor((rows(draws) + 1) :/ (100 \ 20 \ 10))
    cutoffs = J(3, 1, .)
    draws = sort(draws, -1)
    for (j = 1; j <= 3; j++) {
        if (rank[j] > 0) cutoffs[j] = draws[rank[j]]
    }
    return(cutoffs)
}
end
