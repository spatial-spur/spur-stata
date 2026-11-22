*! version 1.1.0  17nov2026 spatial-spur
program define spur, rclass sortpreserve
	version 14

	syntax varlist(numeric min=2) [if] [in], ///
		[ ///
			q(int 15) ///
			nrep(int 100000) ///
			seed(int 42) ///
			latlong ///
			avc(real -1) ///
			uncond ///
			cvs ///
			Replace ///
			BACKend(string) ///
			nthreads(integer -1) ///
		]

	// The combined command needs the backend and result interface provided by scpc-stata.
	capture scpc, capabilities
	if _rc {
		di as error "spur requires a compatible scpc-stata build; install the matching scpc-stata checkout and restart Stata"
		exit 199
	}

	if `q' < 2 {
		di as error "q() must be at least 2"
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

	if `avc' == -1 {
		local avc 0.03
	}
	if missing(`avc') | `avc' <= 0.001 | `avc' >= 0.99 {
		di as error "avc() must lie strictly between 0.001 and 0.99"
		exit 198
	}

	_spur_backend resolve, backend(`"`backend'"') nthreads(`nthreads')
	local req_backend "`r(req_backend)'"
	local nthreads_arg "`r(nthreads)'"
	local plugin_ok = r(plugin_ok)
	local requested_reason "`r(backend_reason)'"

	gettoken depvar indepvars : varlist

	local prefix h_
	local transformation lbmgls
	local ncoef : word count `varlist'
	local k_x_raw : word count `indepvars'

	// Validate destination variable names up front when replace is omitted
	if "`replace'" == "" {
		foreach v of local varlist {
			capture confirm new variable `prefix'`v'
			if _rc {
				di as error "`prefix'`v' already defined or invalid name"
				exit 110
			}
		}
	}

	// Diagnostics, transformation, and both regressions use the same complete cases.
	marksample touse_all
	local use_pipeline = (`plugin_ok' == 1 & "`req_backend'" != "fallback")

	local scpc_options avc(`avc') k(`ncoef') backend(`req_backend') nthreads(`nthreads_arg')
	if "`latlong'" != "" {
		local scpc_options `scpc_options' latlong
	}
	if "`uncond'" != "" {
		local scpc_options `scpc_options' uncond
	}
	if "`cvs'" != "" {
		local scpc_options `scpc_options' cvs
	}

	_spur_backend check_coords, touse(`touse_all') `latlong' quiet
	local slist "`r(slist)'"
	local d = r(d)
	if `q' >= r(n) {
		di as error "q() must be smaller than the number of observations"
		exit 198
	}

	local h_tmp_vars ""
	local n_vars_tr = 1 + `k_x_raw'
	forvalues j = 1/`n_vars_tr' {
		tempvar h_tmp_`j'
		qui gen double `h_tmp_`j'' = .
		local h_tmp_vars "`h_tmp_vars' `h_tmp_`j''"
	}

	tempname diagnostics diagnostic_cvs pipe_tied
	local diag_backend "rust"
	local diag_reason ""
	local tr_reason ""
	local diag_tied 0
	local latlong_flag = ("`latlong'" != "")

	if `use_pipeline' == 1 {
		matrix `diagnostics' = J(4, 3, .)
		matrix `diagnostic_cvs' = J(4, 3, .)
		tempfile emat_file
		set seed `seed'
		_spur_backend write_emat, q(`q') nrep(`nrep') file(`"`emat_file'"')

		local plugin_vars "`slist' `varlist' `h_tmp_vars'"
		local pcmd `"spur_pipeline `d' `k_x_raw' 1 `latlong_flag' `q' `nrep' "`nthreads_arg'" `"`emat_file'"' "`diagnostics'" "`diagnostic_cvs'" "`pipe_tied'""'
		capture noisily _spur_backend plugin_call, vars(`plugin_vars') cond(`touse_all') cmd(`"`pcmd'"')
		local rc = _rc
		capture erase `"`emat_file'"'
		if `rc' != 0 {
			exit `rc'
		}

		local diag_tied = scalar(`pipe_tied')
		if "`req_backend'" == "auto" & `diag_tied' {
			local use_pipeline 0
		}
	}

	if `use_pipeline' == 0 {
		local diag_backend "fallback"
		local diag_reason "`requested_reason'"
		local tr_reason "`requested_reason'"
		if `diag_tied' {
			local diag_reason "tied eigenspace"
			local tr_reason "requested fallback"
		}
		set seed `seed'
		mata: spur_mata_pipeline("`varlist'", "`h_tmp_vars'", "`touse_all'", "`slist'", ///
			`latlong_flag', `q', `nrep', "`diagnostics'", "`diagnostic_cvs'")
	}
	matrix rownames `diagnostics' = i0 i1 i0resid i1resid
	matrix colnames `diagnostics' = teststat p ha_param

	quietly regress `depvar' `indepvars' if `touse_all', robust

	tempname levels_N levels_r2 levels_r2_a
	scalar `levels_N' = e(N)
	scalar `levels_r2' = e(r2)
	scalar `levels_r2_a' = e(r2_a)

	quietly scpc, `scpc_options'

	tempname levels_scpcstats levels_scpccvs
	matrix `levels_scpcstats' = e(scpcstats)
	if "`cvs'" != "" {
		matrix `levels_scpccvs' = e(scpccvs)
	}
	local levels_scpc_be "`e(backend)'"
	local levels_scpc_reason "`e(backend_reason)'"
	local levels_scpc_tied = e(tied)

	local transformed_depvar `prefix'`depvar'
	local transformed_indepvars
	foreach var of local indepvars {
		local transformed_indepvars `transformed_indepvars' `prefix'`var'
	}

	forvalues j = 1/`n_vars_tr' {
		local vname : word `j' of `varlist'
		local hv : word `j' of `h_tmp_vars'
		if "`replace'" != "" {
			capture drop `prefix'`vname'
		}
		rename `hv' `prefix'`vname'
	}

	quietly regress `transformed_depvar' `transformed_indepvars' if `touse_all', robust

	tempname transformed_N transformed_r2 transformed_r2_a
	scalar `transformed_N' = e(N)
	scalar `transformed_r2' = e(r2)
	scalar `transformed_r2_a' = e(r2_a)

	quietly scpc, `scpc_options'

	tempname transformed_scpcstats transformed_scpccvs
	matrix `transformed_scpcstats' = e(scpcstats)
	if "`cvs'" != "" {
		matrix `transformed_scpccvs' = e(scpccvs)
	}
	local tr_scpc_be "`e(backend)'"
	local tr_scpc_reason "`e(backend_reason)'"
	local tr_scpc_tied = e(tied)

	local used_backend "fallback"
	if "`diag_backend'" == "rust" & "`levels_scpc_be'" == "rust" & "`tr_scpc_be'" == "rust" {
		local used_backend "rust"
	}
	local any_tied = max(`diag_tied', `levels_scpc_tied', `tr_scpc_tied')

	display as text "SPUR Diagnostics"
	matlist `diagnostics', border(all) format(%9.4f)

	display as text ""
	display as text "Levels SCPC"
	matlist `levels_scpcstats', border(all)

	display as text ""
	display as text "Transformed SCPC"
	matlist `transformed_scpcstats', border(all)

	local row 0
	foreach diagnostic in i0 i1 i0resid i1resid {
		local ++row
		tempname cv
		matrix `cv' = (`diagnostic_cvs'[`row', 1] \ `diagnostic_cvs'[`row', 2] \ `diagnostic_cvs'[`row', 3])
		return scalar `diagnostic'_teststat = `diagnostics'[`row', 1]
		return scalar `diagnostic'_p = `diagnostics'[`row', 2]
		return scalar `diagnostic'_ha_param = `diagnostics'[`row', 3]
		return matrix `diagnostic'_cv = `cv'
		return local `diagnostic'_backend "`diag_backend'"
		return local `diagnostic'_backend_reason "`diag_reason'"
	}

	return matrix diagnostics = `diagnostics'

	return scalar levels_N = `levels_N'
	return scalar levels_r2 = `levels_r2'
	return scalar levels_r2_a = `levels_r2_a'
	return matrix levels_scpcstats = `levels_scpcstats'

	return scalar transformed_N = `transformed_N'
	return scalar transformed_r2 = `transformed_r2'
	return scalar transformed_r2_a = `transformed_r2_a'
	return matrix transformed_scpcstats = `transformed_scpcstats'

	if "`cvs'" != "" {
		return matrix levels_scpccvs = `levels_scpccvs'
		return matrix transformed_scpccvs = `transformed_scpccvs'
	}

	return scalar tied = `any_tied'
	return local backend "`used_backend'"
	return local diagnostics_backend "`diag_backend'"
	return local diagnostics_backend_reason "`diag_reason'"
	return local transform_backend "`diag_backend'"
	return local transform_backend_reason "`tr_reason'"
	return local scpc_levels_backend "`levels_scpc_be'"
	return local scpc_levels_backend_reason "`levels_scpc_reason'"
	return local scpc_transformed_backend "`tr_scpc_be'"
	return local scpc_transformed_backend_reason "`tr_scpc_reason'"
	return local depvar "`depvar'"
	return local indepvars "`indepvars'"
	return local transformed_depvar "`transformed_depvar'"
	return local transformed_indepvars "`transformed_indepvars'"
	return local prefix "`prefix'"
	return local transformation "`transformation'"
end
