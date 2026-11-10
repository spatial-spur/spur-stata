* Run from the repository root in licensed Stata.
clear all
set more off
adopath ++ .
spur_mata
mata:
sv2_latlongflag = 0
S = (1::20), (1::20):^2
D = sv2_getdistmat(S)
Dcopy = D
B = sv2_demeanmat(D)
assert(max(abs(D - Dcopy)) == 0)
c = sv2_getcbar(0.03, D)
assert(abs(mean(exp(-c * sv2_lvech(D))) - 0.03) < 1e-10)
R = sv2_get_R(sv2_get_sigma_lbm_dm(D), 4)
assert(max(abs(cross(R, R) - I(4))) < 1e-10)
assert(max(abs(colsum(R))) < 1e-8)
sv2_latlongflag = 1
D = sv2_getdistmat((35.123, 45.987 \ -35.123, -134.013))
assert(abs(D[1,2] - 0.5) < 1e-8)
end

capture mata: sv2_get_ha_parm_I0(I(3), I(3), I(3), J(3, 100, 1))
assert _rc == 3498
