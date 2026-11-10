*! version 1.0.0 10nov2026
program define spur_mata
    version 14
    capture mata: assert(sv2_version() == 20261110)
    if _rc {
        quietly findfile spur_mata_source.do
        quietly run "`r(fn)'"
    }
end
