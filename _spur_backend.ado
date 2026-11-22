*! version 1.1.0  17nov2026 spatial-spur
mata mata set matastrict on

program define _spur_backend, rclass
    version 16.0

    gettoken subcmd rest : 0, parse(" ,")
    local rest = ltrim(`"`rest'"')

    if "`subcmd'" == "resolve" {
        local 0 `"`rest'"'
        syntax , [BACKend(string) nthreads(integer -1)]

        local req_backend = lower(trim(`"`backend'"'))
        if "`req_backend'" == "" {
            local req_backend = lower(trim(`"$SPATIAL_BACKEND"'))
        }
        if "`req_backend'" == "" {
                local env_backend : environment SPATIAL_BACKEND
                local req_backend = lower(trim(`"`env_backend'"'))
        }
        if "`req_backend'" == "" {
            local req_backend "auto"
        }
        if "`req_backend'" != "auto" & "`req_backend'" != "rust" & "`req_backend'" != "fallback" {
            di as error "invalid backend() option: `req_backend' (expected auto, rust, or fallback)"
            exit 198
        }
        if `nthreads' != -1 & `nthreads' < 1 {
            di as error "nthreads() must be a positive integer (or -1 for default)"
            exit 198
        }

        local plugin_ok 0
        local backend_reason "requested fallback"
        if "`req_backend'" != "fallback" {
            local plugin_reason "compatible spur_rust.plugin (v4) could not be loaded. Install the matching plugin"
            capture plugin call _spur_rust_plugin, ping_spur_v4
            if _rc == 0 {
                local plugin_ok 1
            }
            else {
                capture plugin call _spur_rust_plugin, ping
                if _rc != 0 {
                    capture program _spur_rust_plugin, plugin using("spur_rust.plugin")
                    if _rc != 0 {
                        capture findfile "spur_rust.plugin"
                        if _rc == 0 {
                            capture program _spur_rust_plugin, plugin using("`r(fn)'")
                        }
                    }
                }
                capture plugin call _spur_rust_plugin, ping_spur_v4
                if _rc == 0 {
                    local plugin_ok 1
                }
                else {
                    capture plugin call _spur_rust_plugin, ping
                    if _rc == 0 {
                        local plugin_reason "loaded spur_rust.plugin is incompatible with the required API (v4). Install the matching plugin and restart Stata"
                    }
                }
            }
            local backend_reason ""
            if `plugin_ok' == 0 {
                local backend_reason "`plugin_reason'"
                if "`req_backend'" == "rust" {
                    di as error "`backend_reason'"
                    exit 199
                }
            }
        }

        return local backend_reason "`backend_reason'"
        return local req_backend "`req_backend'"
        return local nthreads "`nthreads'"
        return scalar plugin_ok = `plugin_ok'
        exit
    }

    if "`subcmd'" == "check_coords" {
        local 0 `"`rest'"'
        syntax , touse(varname) [latlong quiet]

        local slist ""
        local d 0
        capture unab raw_slist : s_*
        if _rc != 0 {
            di as error "variable s_* not found"
            exit 111
        }
        local d : word count `raw_slist'
        forvalues s_idx = 1/`d' {
            local in_list : list posof "s_`s_idx'" in raw_slist
            if `in_list' == 0 {
                disp as text "s_* variables not continuously numbered starting from 1 (s_1, s_2, etc.)"
                exit 999
            }
            capture confirm numeric variable s_`s_idx', exact
            if _rc != 0 {
                disp as text "s_* variables must be numeric (s_1, s_2, etc.)"
                exit 999
            }
            local slist "`slist' s_`s_idx'"
        }
        local slist = ltrim("`slist'")

        foreach svar of local slist {
            qui count if missing(`svar') & `touse'
            if r(N) > 0 {
                disp as text "missing value(s) in s_* variable; aborting"
                exit 999
            }
        }

        qui count if `touse'
        local ns = r(N)
        if `ns' < 5 {
            disp as text "too few locations found in variables s_1, s_2 etc; aborting"
            exit 999
        }

        if "`latlong'" != "" {
            if `d' != 2 {
                di as error "latlong requires exactly two coordinates in s_1 and s_2"
                exit 198
            }
            quietly count if `touse' & (!inrange(s_1, -90, 90) | !inrange(s_2, -180, 180))
            if r(N) > 0 {
                di as error "latlong requires latitude in [-90, 90] and longitude in [-180, 180]"
                exit 198
            }
        }

        if "`quiet'" == "" {
            _spur_backend print_coords, n(`ns') d(`d') `latlong'
        }

        return local slist "`slist'"
        return scalar d = `d'
        return scalar n = `ns'
        exit
    }

    if "`subcmd'" == "print_coords" {
        local 0 `"`rest'"'
        syntax , n(integer) d(integer) [latlong]
        local ns_str = strofreal(`n', "%6.0f")
        local d_str = strofreal(`d', "%3.0f")
        disp as text "found `ns_str' observations and `d_str'-dimensional locations in s_*"
        if "`latlong'" == "" {
            disp as text "using Euclidean norm to compute distance between locations stored in s_*"
        }
        else {
            disp as text "Computing distances on surface of sphere treating s_1 as latitude and s_2 as longitude"
        }
        exit
    }

    if "`subcmd'" == "write_emat" {
        local 0 `"`rest'"'
        syntax , q(integer) nrep(integer) FILE(string)
        capture erase `"`file'"'
        capture noisily {
            mata: _sp_e_E = rnormal(`q', `nrep', 0, 1)
            mata: _sp_e_fh = fopen(`"`file'"', "w")
            mata: _sp_e_C = bufio()
            mata: bufbyteorder(_sp_e_C, 2)
            mata: fbufput(_sp_e_C, _sp_e_fh, "%8z", rowshape(_sp_e_E, 1))
            mata: fclose(_sp_e_fh)
        }
        local rc = _rc
        if `rc' != 0 capture mata: fclose(_sp_e_fh)
        capture mata: mata drop _sp_e_E _sp_e_fh _sp_e_C
        if `rc' != 0 {
            exit `rc'
        }
        exit
    }

    if "`subcmd'" == "plugin_call" {
        local 0 `"`rest'"'
        syntax , vars(string) cond(varname) cmd(string)
        plugin call _spur_rust_plugin `vars' if `cond', `cmd'
        exit
    }

    if "`subcmd'" == "run_test" {
        local 0 `"`rest'"'
        syntax , flavor(string) yvar(string) touse(varname) q(integer) nrep(integer) seed(integer) ///
            [xvars(string) latlong BACKend(string) nthreads(integer -1)]

        if `q' < 2 {
            di as error "q() must be at least 2 for SPUR inference"
            exit 198
        }
        if `nrep' < 1 {
            di as error "nrep() must be at least 1"
            exit 198
        }
        if `seed' < 0 | `seed' >= 4294967296 {
            di as error "seed() must be an integer in [0, 2^32)"
            exit 198
        }

        if "`flavor'" == "i0" {
            local plugin_subcmd "spur_i0"
            local include_int 0
        }
        else if "`flavor'" == "i1" {
            local plugin_subcmd "spur_i1"
            local include_int 0
        }
        else if "`flavor'" == "i0resid" {
            local plugin_subcmd "spur_i0_resid"
            local include_int 1
        }
        else if "`flavor'" == "i1resid" {
            local plugin_subcmd "spur_i1_resid"
            local include_int 1
        }
        else {
            di as error "unknown SPUR test flavor: `flavor'"
            exit 198
        }

        _spur_backend resolve, backend(`"`backend'"') nthreads(`nthreads')
        local req_backend "`r(req_backend)'"
        local nthreads_arg "`r(nthreads)'"
        local plugin_ok = r(plugin_ok)
        local backend_reason "`r(backend_reason)'"

        _spur_backend check_coords, touse(`touse') `latlong' quiet
        local slist "`r(slist)'"
        local d = r(d)
        local n_obs = r(n)

        if `q' >= `n_obs' {
            di as error "q() must be smaller than the number of observations"
            exit 198
        }

        tempname p teststat ha cv tied_scal
        local used_backend ""
        local tied_val 0

        if "`req_backend'" == "fallback" | `plugin_ok' == 0 {
            set seed `seed'
            _spur_backend mata_test, flavor(`flavor') yvar(`yvar') xvars(`xvars') touse(`touse') `latlong' ///
                q(`q') nrep(`nrep') pmat(`p') statscal(`teststat') hascal(`ha') cvmat(`cv')
            local used_backend "fallback"
        }
        else {
            tempfile emat_file
            matrix `cv' = J(3, 1, .)
            set seed `seed'
            _spur_backend write_emat, q(`q') nrep(`nrep') file(`"`emat_file'"')

            local k_x_raw : word count `xvars'
            local latlong_flag = cond("`latlong'" != "", 1, 0)
            local plugin_vars "`slist' `yvar' `xvars'"
            capture noisily plugin call _spur_rust_plugin `plugin_vars' if `touse', ///
                `plugin_subcmd' `d' `k_x_raw' `include_int' `latlong_flag' `q' `nrep' "`nthreads_arg'" ///
                `"`emat_file'"' "`cv'" "`teststat'" "`p'" "`ha'" "`tied_scal'"
            local rc = _rc
            capture erase `"`emat_file'"'
            if `rc' != 0 {
                exit `rc'
            }

            if "`req_backend'" == "auto" & scalar(`tied_scal') == 1 {
                set seed `seed'
                _spur_backend mata_test, flavor(`flavor') yvar(`yvar') xvars(`xvars') touse(`touse') `latlong' ///
                    q(`q') nrep(`nrep') pmat(`p') statscal(`teststat') hascal(`ha') cvmat(`cv')
                local used_backend "fallback"
                local tied_val 1
                local backend_reason "tied eigenspace"
            }
            else {
                _spur_backend print_coords, n(`n_obs') d(`d') `latlong'
                local used_backend "rust"
                local tied_val = scalar(`tied_scal')
            }
        }

        return scalar p = scalar(`p')
        return scalar teststat = scalar(`teststat')
        return scalar ha_param = scalar(`ha')
        return matrix cv = `cv'
        return scalar tied = `tied_val'
        return local backend "`used_backend'"
        return local backend_reason "`backend_reason'"
        exit
    }

    if "`subcmd'" == "mata_test" {
        local 0 `"`rest'"'
        syntax , flavor(string) yvar(string) touse(varname) q(integer) nrep(integer) ///
            pmat(name) statscal(name) hascal(name) cvmat(name) [xvars(string) latlong]

        mata: spur_mata_test("`flavor'", "`yvar'", "`xvars'", "`touse'", "`latlong'", ///
            `q', `nrep', "`pmat'", "`statscal'", "`hascal'", "`cvmat'")
        exit
    }

    if "`subcmd'" == "mata_halflife" {
        local 0 `"`rest'"'
        syntax , yvar(string) touse(varname) q(integer) nrep(integer) level(real) ///
            cil(name) ciu(name) maxdist(name) [latlong]

        capture noisily {
            mata: _sp_t_e = rnormal(`q', `nrep', 0, 1)
            mata: _sp_t_y = st_data(., "`yvar'", "`touse'")
            mata: _sp_t_s = get_s_matrix("`touse'", "`latlong'")
            mata: _sp_t_diam = .
            mata: _sp_t_d = getdistmat_normalized(_sp_t_s, _sp_t_diam)
            mata: _sp_t_r = spatial_persistence(_sp_t_y, _sp_t_d, _sp_t_e, `level' / 100)
            mata: st_numscalar("`cil'", _sp_t_r.hl_ci_lower)
            mata: st_numscalar("`ciu'", _sp_t_r.hl_ci_upper)
            if "`latlong'" != "" {
                mata: st_numscalar("`maxdist'", _sp_t_diam * 2 * pi() * 6371000.009)
            }
            else {
                mata: st_numscalar("`maxdist'", _sp_t_diam)
            }
        }
        local rc = _rc
        capture mata: mata drop _sp_t_e _sp_t_y _sp_t_s _sp_t_d _sp_t_diam _sp_t_r
        if `rc' != 0 {
            exit `rc'
        }
        exit
    }

    if "`subcmd'" == "mata_transform" {
        local 0 `"`rest'"'
        syntax , invars(string) outvars(string) touse(varname) ///
            transformation(string) radius(string) [clustvar(string) latlong]

        capture noisily {
            mata: _sp_t_s = get_s_matrix("`touse'", "`latlong'")
            mata: _sp_t_y = st_data(., "`invars'", "`touse'")
            if "`transformation'" == "lbmgls" {
                mata: _sp_t_hy = spur_lbm_gls(_sp_t_s, _sp_t_y)
            }
            else {
                if "`transformation'" == "cluster" {
                    mata: _sp_t_cl = get_cluster_matrix("`clustvar'", "`touse'")
                }
                else {
                    mata: _sp_t_cl = 0
                }
                mata: _sp_t_hy = make_transform(_sp_t_s, "`transformation'", strtoreal("`radius'"), _sp_t_cl) * _sp_t_y
            }
            mata: st_store(., tokens("`outvars'"), "`touse'", _sp_t_hy)
        }
        local rc = _rc
        capture mata: mata drop _sp_t_s _sp_t_cl _sp_t_y _sp_t_hy
        if `rc' != 0 {
            exit `rc'
        }
        exit
    }

    di as error "unknown _spur_backend subcommand: `subcmd'"
    exit 198
end
