*! version 1.0.0  21jan2025 pdavidboll
version 16.0
mata

real scalar getpow_qf(real matrix om0, real matrix om1, real matrix e) {
    real matrix lower0, lower1, ho, ha, ya_o, yo_a
    real vector qe, qa_o, qo_a, lr_o, lr_a
    real scalar cv, pow

    lower0 = cholesky(om0)
    lower1 = cholesky(om1)
    ho = lusolve(lower1, lower0)
    ha = lusolve(lower0, lower1)

    // Step 4: Compute quadratic forms
    qe = colsum(e:^2)    // Sum of squares of columns of e
    ya_o = ho * e
    yo_a = ha * e
    qa_o = colsum(ya_o:^2) // Sum of squares of columns of ya_o
    qo_a = colsum(yo_a:^2) // Sum of squares of columns of yo_a
  
    // Step 5: Likelihood ratios
    lr_o = qe :/ qa_o
    lr_a = qo_a :/ qe
  
    // Step 6: Compute 95th percentile critical value and calculate power
    cv = spur_quantile(lr_o', 0.95) // 95th percentile of lr_o
    pow = mean(lr_a' :> cv)
  
    return(pow)
}

end