version 14
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
end
