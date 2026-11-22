version 16.0
mata
real scalar spur_shape_equal(real matrix a, real matrix b)
{
    real scalar ta, tb
    ta = trace(a)
    tb = trace(b)
    if (ta <= 0 | tb <= 0 | ta >= . | tb >= .) {
        errprintf("SPUR covariance has non-positive trace\n")
        exit(198)
    }
    return(max(abs(a :/ ta - b :/ tb)) <= 1e-12)
}
end
