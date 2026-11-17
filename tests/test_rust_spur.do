*! Comprehensive Stata Rust & Mata Fallback parity and hygiene suite for spur-stata.
clear all
mata mata clear
set more off
version 14

local repo_root "/Users/pauldavidboll/Github/spatial-spur"
adopath ++ "`repo_root'/scpc-stata/src"
adopath ++ "`repo_root'/spur-stata"

// Preload scpc.ado once before defining test Mata functions, since scpc.ado
// executes `mata mata clear` on initial load.
capture scpc

capture program drop assert_mat_close
program define assert_mat_close
    args m1 m2 tol
    assert rowsof(`m1') == rowsof(`m2')
    assert colsof(`m1') == colsof(`m2')
    local nr = rowsof(`m1')
    local nc = colsof(`m1')
    forvalues i = 1/`nr' {
        forvalues j = 1/`nc' {
            local v1 = `m1'[`i', `j']
            local v2 = `m2'[`i', `j']
            if missing(`v1') | missing(`v2') {
                assert missing(`v1') & missing(`v2')
            }
            else {
                assert abs(`v1' - `v2') < `tol'
            }
        }
    }
end

capture program drop assert_var_close
program define assert_var_close
    args v1 v2 tol
    tempvar miss_mismatch diff
    qui gen byte `miss_mismatch' = (missing(`v1') != missing(`v2'))
    qui count if `miss_mismatch'
    assert r(N) == 0
    qui gen double `diff' = abs(`v1' - `v2') if !missing(`v1')
    qui summarize `diff', meanonly
    if r(N) > 0 {
        assert r(max) < `tol'
    }
end

mata:
mata set matastrict on

real colvector _test_canon_signs(real matrix R)
{
    real scalar j, i, best_i, best_abs, v
    real colvector s
    s = J(cols(R), 1, 1)
    for (j = 1; j <= cols(R); j++) {
        best_i = 1
        best_abs = abs(R[1, j])
        for (i = 2; i <= rows(R); i++) {
            v = abs(R[i, j])
            if (v > best_abs) {
                best_abs = v
                best_i = i
            }
        }
        if (R[best_i, j] < 0) s[j] = -1
    }
    return(s)
}

void _test_verify_canon_test(
    string scalar flavor,
    string scalar yname,
    string scalar xnames,
    string scalar touse,
    string scalar latlong,
    real scalar q,
    real scalar nrep,
    real scalar rust_p,
    real scalar rust_stat,
    real scalar rust_ha,
    string scalar rust_cv_mat,
    real scalar cv_tol
)
{
    real matrix y, X, s, distmat, sig, R, emat, emat_al, cv_r, M
    real colvector sg
    real scalar c_bm, cv_diff
    struct struct_spatial_i0_test scalar t0
    struct struct_spatial_i1_test scalar t1
    struct struct_spatial_i0_test scalar t0r
    struct struct_spatial_i1_test_residual scalar t1r

    emat = rnormal(q, nrep, 0, 1)
    y = st_data(., yname, touse)
    s = get_s_matrix(touse, latlong)
    distmat = getdistmat_normalized(s)

    if (flavor == "i0" | flavor == "i1") {
        sig = get_sigma_lbm_dm(distmat)
        R = get_R(sig, q)
        sg = _test_canon_signs(R)
        emat_al = sg :* emat
        if (flavor == "i0") {
            t0 = spatial_i0_test(y, distmat, emat_al)
            assert(abs(t0.pvalue[1] - rust_p) < 1e-10)
            assert(abs(t0.LR[1] - rust_stat) < 1e-10 * max((1, abs(rust_stat))))
            assert(abs(t0.ha_parm - rust_ha) < 1e-10 * max((1, abs(rust_ha))))
            cv_r = st_matrix(rust_cv_mat)
            cv_diff = max(abs(t0.cvalue - cv_r) :/ (1 :+ abs(cv_r)))
            assert(cv_diff < cv_tol)
        }
        else {
            t1 = spatial_i1_test(y, distmat, emat_al)
            assert(abs(t1.pvalue[1] - rust_p) < 1e-10)
            assert(abs(t1.LR[1] - rust_stat) < 1e-10 * max((1, abs(rust_stat))))
            assert(abs(t1.ha_parm - rust_ha) < 1e-10 * max((1, abs(rust_ha))))
            cv_r = st_matrix(rust_cv_mat)
            cv_diff = max(abs(t1.cv_vec - cv_r) :/ (1 :+ abs(cv_r)))
            assert(cv_diff < cv_tol)
        }
    }
    else {
        if (strtrim(xnames) != "") {
            X = st_data(., xnames, touse)
            X = J(rows(X), 1, 1), X
        }
        else {
            X = J(rows(y), 1, 1)
        }
        M = I(rows(y)) - X * invsym(X' * X) * X'
        c_bm = getcbar(0.999, distmat)
        sig = get_sigma_residual(distmat, c_bm, M)
        R = get_R(sig, q)
        sg = _test_canon_signs(R)
        emat_al = sg :* emat
        if (flavor == "i0resid") {
            t0r = spatial_i0_test_residual(y, X, distmat, emat_al)
            assert(abs(t0r.pvalue[1] - rust_p) < 1e-10)
            assert(abs(t0r.LR[1] - rust_stat) < 1e-10 * max((1, abs(rust_stat))))
            assert(abs(t0r.ha_parm - rust_ha) < 1e-10 * max((1, abs(rust_ha))))
            cv_r = st_matrix(rust_cv_mat)
            cv_diff = max(abs(t0r.cvalue - cv_r) :/ (1 :+ abs(cv_r)))
            assert(cv_diff < cv_tol)
        }
        else {
            t1r = spatial_i1_test_residual(y, X, distmat, emat_al)
            assert(abs(t1r.pvalue[1] - rust_p) < 1e-10)
            assert(abs(t1r.LR[1] - rust_stat) < 1e-10 * max((1, abs(rust_stat))))
            assert(abs(t1r.ha_parm - rust_ha) < 1e-10 * max((1, abs(rust_ha))))
            cv_r = st_matrix(rust_cv_mat)
            cv_diff = max(abs(t1r.cv_vec - cv_r) :/ (1 :+ abs(cv_r)))
            assert(cv_diff < cv_tol)
        }
    }
}

void _test_verify_canon_halflife(
    string scalar yname,
    string scalar touse,
    string scalar latlong,
    real scalar q,
    real scalar nrep,
    real scalar level,
    real scalar normdist,
    real scalar rust_cil,
    real scalar rust_ciu,
    real scalar rust_md
)
{
    real matrix y, s, distmat, sig, R, emat, emat_al
    real colvector sg
    real scalar c_bm, md, cil, ciu
    struct struct_spatial_persistence scalar ci

    emat = rnormal(q, nrep, 0, 1)
    y = st_data(., yname, touse)
    s = get_s_matrix(touse, latlong)
    distmat = getdistmat_normalized(s)
    c_bm = getcbar(0.999, distmat)
    sig = get_sigma_dm(distmat, c_bm)
    R = get_R(sig, q)
    sg = _test_canon_signs(R)
    emat_al = sg :* emat

    ci = spatial_persistence(y, distmat, emat_al, level / 100)
    md = max(getdistmat(s)) * 3.14159265359 * 6371000.009 * 2
    assert(abs(md - rust_md) < 1e-6)
    cil = ci.hl_ci_lower
    ciu = ci.hl_ci_upper
    if (ciu >= 100) ciu = .
    if (!normdist) {
        cil = cil * md
        ciu = ciu * md
    }
    assert(abs(cil - rust_cil) < 1e-10 * max((1, abs(cil))))
    if (missing(ciu) | missing(rust_ciu)) {
        assert(missing(ciu) & missing(rust_ciu))
    }
    else {
        assert(abs(ciu - rust_ciu) < 1e-10 * max((1, abs(ciu))))
    }
}
end

// ---------------------------------------------------------------------------
// 1. Geodesic (latlong) parity on bench_200.csv (s_1=lat, s_2=lon)
// ---------------------------------------------------------------------------
import delimited "`repo_root'/verification/data/bench_200.csv", clear asdouble
drop s_1_clustered s_2_clustered
rename s_1 _tmp_lon
rename s_2 s_1
rename _tmp_lon s_2

// Preserve user variables, sort order, and Mata globals across all commands
gen long obs_order = _n
gen double user_keep_var = sin(_n)
mata: emat = (999)
mata: y = (888)
mata: H = (777)

// Shuffle observation order to exercise sortpreserve
set seed 123
gen double _rand_sort = runiform()
sort _rand_sort
drop _rand_sort
gen long shuffled_order = _n
gen byte _all_one = 1

// 1a. All four diagnostics: exact canonical-sign parity & public fallback execution
foreach flavor in i0 i1 i0resid i1resid {
    local vlist "y"
    local xrest ""
    if "`flavor'" == "i0resid" | "`flavor'" == "i1resid" {
        local vlist "y x1 x2"
        local xrest "x1 x2"
    }

    spurtest `flavor' `vlist', q(12) nrep(401) seed(42) latlong backend(rust)
    assert "`r(backend)'" == "rust"
    assert r(tied) == 0
    local r_p = r(p)
    local r_stat = r(teststat)
    local r_ha = r(ha_param)
    matrix r_cv = r(cv)

    set seed 42
    // cv uses Hyndman-Fan type 7 in Rust vs mm_quantile def 2 in Mata (O(1/nrep) tail spacing)
    mata: _test_verify_canon_test("`flavor'", "y", "`xrest'", "_all_one", "latlong", 12, 401, `r_p', `r_stat', `r_ha', "r_cv", 0.05)

    // Public fallback also runs cleanly without altering user Mata globals
    spurtest `flavor' `vlist', q(12) nrep(401) seed(42) latlong backend(fallback)
    assert "`r(backend)'" == "fallback"
}

// Verify tight cv convergence between def 7 (Rust) and def 2 (Mata) at large nrep=20000
spurtest i0 y, q(15) nrep(20000) seed(42) latlong backend(rust)
local r_p = r(p)
local r_stat = r(teststat)
local r_ha = r(ha_param)
matrix r_cv = r(cv)
set seed 42
mata: _test_verify_canon_test("i0", "y", "", "_all_one", "latlong", 15, 20000, `r_p', `r_stat', `r_ha', "r_cv", 1e-3)

// Also test public fallback == rust at default nrep=100000 for i0 (where eigensigns don't shift ha bin)
spurtest i0 y, q(15) nrep(100000) seed(42) latlong backend(rust)
local r_stat = r(teststat)
local r_ha = r(ha_param)
spurtest i0 y, q(15) nrep(100000) seed(42) latlong backend(fallback)
assert abs(r(teststat) - `r_stat') < 1e-10
assert abs(r(ha_param) - `r_ha') < 1e-10

// Residual test with single intercept regressor (no covariates)
foreach flavor in i0resid i1resid {
    spurtest `flavor' y, q(10) nrep(301) seed(7) latlong backend(rust)
    assert "`r(backend)'" == "rust"
    local r_p = r(p)
    local r_stat = r(teststat)
    local r_ha = r(ha_param)
    matrix r_cv = r(cv)
    set seed 7
    mata: _test_verify_canon_test("`flavor'", "y", "", "_all_one", "latlong", 10, 301, `r_p', `r_stat', `r_ha', "r_cv", 0.05)
}

// 1b. If/in subsetting parity
spurtest i0 y if _n <= 110 in 10/130, q(8) nrep(301) seed(11) latlong backend(rust)
assert "`r(backend)'" == "rust"
local r_p = r(p)
local r_stat = r(teststat)
local r_ha = r(ha_param)
matrix r_cv = r(cv)
gen byte _sub = (_n >= 10 & _n <= 130 & _n <= 110)
set seed 11
mata: _test_verify_canon_test("i0", "y", "", "_sub", "latlong", 8, 301, `r_p', `r_stat', `r_ha', "r_cv", 0.05)
drop _sub

// 1c. Half-life CI exact canonical parity + public fallback execution (latlong, meters & normdist, 95 & 90)
spurhalflife y, q(8) nrep(301) seed(42) level(95) latlong backend(rust)
assert "`r(backend)'" == "rust"
local r_cil = r(ci_l)
local r_ciu = r(ci_u)
local r_md = r(max_dist)
set seed 42
mata: _test_verify_canon_halflife("y", "_all_one", "latlong", 8, 301, 95, 0, `r_cil', `r_ciu', `r_md')
spurhalflife y, q(8) nrep(301) seed(42) level(95) latlong backend(fallback)
assert "`r(backend)'" == "fallback"

spurhalflife y, q(8) nrep(301) seed(42) level(90) latlong normdist backend(rust)
assert "`r(backend)'" == "rust"
local r_cil = r(ci_l)
local r_ciu = r(ci_u)
local r_md = r(max_dist)
set seed 42
mata: _test_verify_canon_halflife("y", "_all_one", "latlong", 8, 301, 90, 1, `r_cil', `r_ciu', `r_md')

// ---------------------------------------------------------------------------
// 2. Euclidean parity on bench_200.csv (i0, i1, i0resid, i1resid, halflife)
// ---------------------------------------------------------------------------
foreach flavor in i0 i1 i0resid i1resid {
    local vlist "y"
    local xrest ""
    if "`flavor'" == "i0resid" | "`flavor'" == "i1resid" {
        local vlist "y x1 x2"
        local xrest "x1 x2"
    }
    spurtest `flavor' `vlist', q(10) nrep(301) seed(19) backend(rust)
    assert "`r(backend)'" == "rust"
    local r_p = r(p)
    local r_stat = r(teststat)
    local r_ha = r(ha_param)
    matrix r_cv = r(cv)
    set seed 19
    mata: _test_verify_canon_test("`flavor'", "y", "`xrest'", "_all_one", "", 10, 301, `r_p', `r_stat', `r_ha', "r_cv", 0.05)
}

spurhalflife y, q(8) nrep(301) seed(19) level(95) normdist backend(rust)
local r_cil = r(ci_l)
local r_ciu = r(ci_u)
local r_md = r(max_dist)
set seed 19
mata: _test_verify_canon_halflife("y", "_all_one", "", 8, 301, 95, 1, `r_cil', `r_ciu', `r_md')

drop _all_one

// ---------------------------------------------------------------------------
// 3. Spatial transformations parity (lbmgls, nn, iso with isolated pts, cluster)
// ---------------------------------------------------------------------------
// Note: Mata getdistmat uses 12-digit 3.14159265359 vs IEEE pi (~8.5e-13 rel),
// which sigma_lbm^{-1/2} amplifies up to ~1e-7 on dense latlong points.
spurtransform y x1 x2, prefix(rl_) transformation(lbmgls) latlong backend(rust)
assert "`r(backend)'" == "rust"
spurtransform y x1 x2, prefix(fl_) transformation(lbmgls) latlong backend(fallback)
assert "`r(backend)'" == "fallback"
foreach v in y x1 x2 {
    assert_var_close rl_`v' fl_`v' 1e-6
}
drop rl_* fl_*

spurtransform y x1, prefix(rn_) transformation(nn) latlong backend(rust)
assert "`r(backend)'" == "rust"
spurtransform y x1, prefix(fn_) transformation(nn) latlong backend(fallback)
assert "`r(backend)'" == "fallback"
assert_var_close rn_y fn_y 1e-10
assert_var_close rn_x1 fn_x1 1e-10
drop rn_* fn_*

// Small radius iso (isolated obs returning 0)
spurtransform y x1, prefix(ri_) transformation(iso) radius(25000) latlong backend(rust)
assert "`r(backend)'" == "rust"
spurtransform y x1, prefix(fi_) transformation(iso) radius(25000) latlong backend(fallback)
assert_var_close ri_y fi_y 1e-10
assert_var_close ri_x1 fi_x1 1e-10
drop ri_* fi_*

// Cluster with string group and subset skipping a cluster
gen str10 grp_str = cond(mod(obs_order, 5) == 0, "G5", cond(mod(obs_order, 3) == 0, "G3", "G1"))
spurtransform y x1 if grp_str != "G3", prefix(rc_) transformation(cluster) clustvar(grp_str) backend(rust)
assert "`r(backend)'" == "rust"
spurtransform y x1 if grp_str != "G3", prefix(fc_) transformation(cluster) clustvar(grp_str) backend(fallback)
assert_var_close rc_y fc_y 1e-10
assert_var_close rc_x1 fc_x1 1e-10
drop rc_* fc_* grp_str

// Separately + missing values + replace + collision error safety
preserve
replace x1 = . in 1/12
spurtransform y x1, prefix(rs_) transformation(lbmgls) latlong separately backend(rust)
spurtransform y x1, prefix(fs_) transformation(lbmgls) latlong separately backend(fallback)
assert_var_close rs_y fs_y 1e-6
assert_var_close rs_x1 fs_x1 1e-6

gen double saved_rs_y = rs_y
capture spurtransform y x1, prefix(rs_) transformation(lbmgls) latlong backend(rust)
assert _rc == 110
assert_var_close rs_y saved_rs_y 1e-15

capture spurtransform y in 1/3, prefix(rs_) transformation(lbmgls) latlong replace backend(rust)
assert _rc != 0
assert_var_close rs_y saved_rs_y 1e-15
restore

// Verify observation order and Mata globals remained untouched
assert shuffled_order == _n
assert abs(user_keep_var - sin(obs_order)) < 1e-15
mata: assert(emat == 999)
mata: assert(y == 888)
mata: assert(H == 777)

// ---------------------------------------------------------------------------
// 4. Full SPUR pipeline thread invariance, levels/transformed SCPC parity & fallback
// ---------------------------------------------------------------------------
spur y x1 x2, q(10) nrep(301) seed(42) latlong cvs backend(rust) nthreads(1)
assert "`r(backend)'" == "rust"
assert r(tied) == 0
matrix r1_diag = r(diagnostics)
matrix r1_i0_cv = r(i0_cv)
matrix r1_i1_cv = r(i1_cv)
matrix r1_i0r_cv = r(i0resid_cv)
matrix r1_i1r_cv = r(i1resid_cv)
matrix r1_lev = r(levels_scpcstats)
matrix r1_tr = r(transformed_scpcstats)
matrix r1_lev_cvs = r(levels_scpccvs)
matrix r1_tr_cvs = r(transformed_scpccvs)
matrix r1_e_stats = e(scpcstats)
assert "`e(backend)'" == "rust"
rename (h_y h_x1 h_x2) (p1_y p1_x1 p1_x2)

// Verify each row of r(diagnostics) matches canonical-sign Mata diagnostics to 1e-10
gen byte _all_one = 1
set seed 42
local r1_p = r1_diag[1, 2]
local r1_st = r1_diag[1, 1]
local r1_ha = r1_diag[1, 3]
mata: _test_verify_canon_test("i0", "y", "", "_all_one", "latlong", 10, 301, `r1_p', `r1_st', `r1_ha', "r1_i0_cv", 0.05)
set seed 42
local r1_p = r1_diag[2, 2]
local r1_st = r1_diag[2, 1]
local r1_ha = r1_diag[2, 3]
mata: _test_verify_canon_test("i1", "y", "", "_all_one", "latlong", 10, 301, `r1_p', `r1_st', `r1_ha', "r1_i1_cv", 0.05)
set seed 42
local r1_p = r1_diag[3, 2]
local r1_st = r1_diag[3, 1]
local r1_ha = r1_diag[3, 3]
mata: _test_verify_canon_test("i0resid", "y", "x1 x2", "_all_one", "latlong", 10, 301, `r1_p', `r1_st', `r1_ha', "r1_i0r_cv", 0.05)
set seed 42
local r1_p = r1_diag[4, 2]
local r1_st = r1_diag[4, 1]
local r1_ha = r1_diag[4, 3]
mata: _test_verify_canon_test("i1resid", "y", "x1 x2", "_all_one", "latlong", 10, 301, `r1_p', `r1_st', `r1_ha', "r1_i1r_cv", 0.05)
drop _all_one

// Thread invariance nthreads(4)
spur y x1 x2, q(10) nrep(301) seed(42) latlong cvs backend(rust) nthreads(4)
assert "`r(backend)'" == "rust"
matrix r4_diag = r(diagnostics)
matrix r4_tr = r(transformed_scpcstats)
assert_mat_close r1_diag r4_diag 1e-12
assert_mat_close r1_tr r4_tr 1e-12
assert_var_close h_y p1_y 1e-12
assert_var_close h_x1 p1_x1 1e-12
assert_var_close h_x2 p1_x2 1e-12

// Verify levels_scpcstats, transformed_scpcstats, h_*, and replace against backend(fallback)
spur y x1 x2, q(10) nrep(301) seed(42) latlong cvs backend(fallback) replace
assert "`r(backend)'" == "fallback"
matrix fb_lev = r(levels_scpcstats)
matrix fb_tr = r(transformed_scpcstats)
matrix fb_lev_cvs = r(levels_scpccvs)
matrix fb_tr_cvs = r(transformed_scpccvs)
assert_mat_close r1_lev fb_lev 1e-10
assert_mat_close r1_tr fb_tr 1e-6
assert_mat_close r1_lev_cvs fb_lev_cvs 1e-10
assert_mat_close r1_tr_cvs fb_tr_cvs 1e-6
assert_mat_close r1_e_stats fb_tr 1e-6
assert_var_close h_y p1_y 1e-6
assert_var_close h_x1 p1_x1 1e-6
assert_var_close h_x2 p1_x2 1e-6
drop p1_*

// Without cvs: ensure levels_scpccvs / transformed_scpccvs absent
spur y x1 x2, q(10) nrep(300) seed(42) latlong backend(rust) replace
capture confirm matrix r(levels_scpccvs)
assert _rc != 0
capture confirm matrix r(transformed_scpccvs)
assert _rc != 0

// Survive user `mata: mata clear` cleanly without losing plugin/emat capability
mata: mata clear
mata: emat = (999)
mata: y = (888)
mata: H = (777)
spurtest i0 y, q(10) nrep(301) seed(42) latlong backend(rust)
assert "`r(backend)'" == "rust"
mata: assert(emat == 999)

// Partial missing in x1 -> different touse between i0(y) and i0resid(y x1)
preserve
replace x1 = . in 1/10
spur y x1 x2, q(10) nrep(250) seed(42) latlong backend(rust) replace
assert "`r(backend)'" == "rust"
matrix pm_tr_r = r(transformed_scpcstats)
rename (h_y h_x1 h_x2) (pm_y pm_x1 pm_x2)
spur y x1 x2, q(10) nrep(250) seed(42) latlong backend(fallback) replace
assert "`r(backend)'" == "fallback"
matrix pm_tr_f = r(transformed_scpcstats)
assert_mat_close pm_tr_r pm_tr_f 1e-6
assert_var_close h_y pm_y 1e-6
assert_var_close h_x1 pm_x1 1e-6
assert_var_close h_x2 pm_x2 1e-6
restore

// $SPATIAL_BACKEND env/global support
global SPATIAL_BACKEND "fallback"
spurtest i0 y, q(10) nrep(200) seed(42) latlong
assert "`r(backend)'" == "fallback"
global SPATIAL_BACKEND ""

// Check shuffled observation order preserved
assert shuffled_order == _n
assert abs(user_keep_var - sin(obs_order)) < 1e-15

// ---------------------------------------------------------------------------
// 5. Tied eigenspaces (regular polygon) auto-fallback vs explicit rust
// ---------------------------------------------------------------------------
clear
set obs 20
gen double theta = 2 * _pi * (_n - 1) / _N
gen double s_1 = cos(theta)
gen double s_2 = sin(theta)
gen double x1 = sin(2 * theta) + 0.1 * cos(3 * theta)
gen double y = cos(theta) + 0.5 * x1

spurtest i0 y, q(15) nrep(200) seed(42) backend(auto)
assert r(tied) == 1
assert "`r(backend)'" == "fallback"

spurtest i0 y, q(15) nrep(200) seed(42) backend(rust)
assert r(tied) == 1
assert "`r(backend)'" == "rust"

spur y x1, q(15) nrep(200) seed(42) backend(auto)
assert r(tied) == 1
assert "`r(backend)'" == "fallback"
rename (h_y h_x1) (auto_hy auto_hx1)
spur y x1, q(15) nrep(200) seed(42) backend(fallback) replace
assert_var_close auto_hy h_y 1e-12
assert_var_close auto_hx1 h_x1 1e-12

// ---------------------------------------------------------------------------
// 6. NN ties & co-located points in lbmgls (rank < n - 1)
// ---------------------------------------------------------------------------
clear
set obs 20
gen double s_1 = mod(_n - 1, 5)
gen double s_2 = floor((_n - 1) / 5)
replace s_1 = s_1[1] in 2
replace s_2 = s_2[1] in 2
gen double y = sin(_n * 0.7)

spurtransform y, prefix(rl_) transformation(lbmgls) backend(rust)
spurtransform y, prefix(fl_) transformation(lbmgls) backend(fallback)
assert_var_close rl_y fl_y 1e-10

spurtransform y, prefix(rn_) transformation(nn) backend(rust)
spurtransform y, prefix(fn_) transformation(nn) backend(fallback)
assert_var_close rn_y fn_y 1e-10

// ---------------------------------------------------------------------------
// 7. Error handling on invalid backend / nthreads / options
// ---------------------------------------------------------------------------
capture spurtest i0 y, backend(bogus)
assert _rc == 198
capture spurtest i0 y, nthreads(0)
assert _rc == 198
capture spurhalflife y, level(100)
assert _rc == 198
capture spurtransform y, prefix(t_) transformation(nn) radius(1)
assert _rc == 198

// Final check on user Mata globals
mata: assert(emat == 999)
mata: assert(y == 888)
mata: assert(H == 777)

display "SPUR Stata Rust tests passed."
