* Run in a fresh Stata session with paths to a legacy plugin and this checkout's package directory.
clear all
set more off
version 14
args oldplugin package_dir scpc_dir
assert `"`oldplugin'"' != "" & `"`package_dir'"' != ""
adopath ++ "`package_dir'"
if `"`scpc_dir'"' != "" {
    adopath ++ "`scpc_dir'"
}
else {
    capture adopath ++ "scpc-stata/src"
}
mata: mata mlib index
program _spur_rust_plugin, plugin using(`"`oldplugin'"')
plugin call _spur_rust_plugin, ping_spur_v1
capture plugin call _spur_rust_plugin, ping_spur_v3
assert _rc != 0
set seed 42
set obs 32
generate double s_1 = runiform()
generate double s_2 = runiform()
generate double x1 = rnormal()
generate double x2 = rnormal()
generate double y = x1 + rnormal()
spurtest i0resid y x1 x2, q(8) nrep(401) seed(42) backend(auto)
assert "`r(backend)'" == "fallback"
assert strpos("`r(backend_reason)'", "(v3)") > 0
assert strpos("`r(backend_reason)'", "loaded spur_rust.plugin is incompatible") > 0
local p = r(p)
local stat = r(teststat)
spurtest i0resid y x1 x2, q(8) nrep(401) seed(42) backend(fallback)
assert abs(r(p) - `p') < 1e-12
assert abs(r(teststat) - `stat') < 1e-10
capture spurtest i0resid y x1 x2, q(8) nrep(401) seed(42) backend(rust)
assert _rc == 199

* Combined spur pipeline with a legacy plugin falls back across all SPUR stages cleanly.
spur y x1 x2, q(8) nrep(401) seed(42) backend(auto)
assert "`r(backend)'" == "fallback"
assert "`r(diagnostics_backend)'" == "fallback"
assert "`r(transform_backend)'" == "fallback"
assert strpos("`r(diagnostics_backend_reason)'", "(v3)") > 0
assert strpos("`r(transform_backend_reason)'", "(v3)") > 0
drop h_y h_x1 h_x2
capture spur y x1 x2, q(8) nrep(401) seed(42) backend(rust)
assert _rc == 199

display "SPUR STALE PLUGIN TEST PASSED"
