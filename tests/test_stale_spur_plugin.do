* Run in a fresh Stata session with paths to a v1-only plugin and isolated package directory.
clear all
set more off
version 14
args oldplugin package_dir
assert `"`oldplugin'"' != ""
local repo_root "/Users/pauldavidboll/Github/spatial-spur"
adopath ++ "`repo_root'/scpc-stata/src"
adopath ++ "`package_dir'"
program _spur_rust_plugin, plugin using(`"`oldplugin'"')
plugin call _spur_rust_plugin, ping_spur_v1
capture plugin call _spur_rust_plugin, ping_spur_v2
assert _rc != 0
import delimited "`repo_root'/verification/data/bench_200.csv", clear asdouble
drop s_1_clustered s_2_clustered
spurtest i0resid y x1 x2, q(8) nrep(401) seed(42) backend(auto)
assert "`r(backend)'" == "fallback"
local p = r(p)
local stat = r(teststat)
spurtest i0resid y x1 x2, q(8) nrep(401) seed(42) backend(fallback)
assert abs(r(p) - `p') < 1e-12
assert abs(r(teststat) - `stat') < 1e-10
capture spurtest i0resid y x1 x2, q(8) nrep(401) seed(42) backend(rust)
assert _rc == 199
display "SPUR STALE PLUGIN TEST PASSED"
