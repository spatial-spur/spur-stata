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

	// check scpc installed
	cap which scpc
	if _rc {
		display as error in smcl `"Please install package {it:scpc} in order to use this command;"' _newline ///
        `"which you can do so by clicking this link: {browse "https://github.com/ukmueller/SCPC/":SCPC}"'
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
	if `avc' < 0.001 | `avc' > 0.99 {
		di as error "avc() must be between 0.001 and 0.99"
		exit 198
	}

	_spur_backend resolve, backend(`"`backend'"') nthreads(`nthreads')
	local req_backend "`r(req_backend)'"
	local nthreads_arg "`r(nthreads)'"
	local plugin_ok = r(plugin_ok)

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

	marksample touse_all
	marksample touse_y, novarlist
	qui replace `touse_y' = 0 if missing(`depvar')
	qui count if `touse_y' != `touse_all'
	local same_sample = (r(N) == 0)

	local use_pipeline = (`plugin_ok' == 1 & "`req_backend'" != "fallback" & `same_sample' == 1)

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

	tempname i0_teststat i0_p i0_ha_param i0_cv
	tempname i1_teststat i1_p i1_ha_param i1_cv
	tempname i0resid_teststat i0resid_p i0resid_ha_param i0resid_cv
	tempname i1resid_teststat i1resid_p i1resid_ha_param i1resid_cv
	tempname diagnostics

	local diag_backend ""
	local diag_tied 0

	if `use_pipeline' == 1 {
		_spur_backend check_coords, touse(`touse_all') `latlong' quiet
		local slist "`r(slist)'"
		local d = r(d)

		local h_tmp_vars ""
		local n_vars_tr = 1 + `k_x_raw'
		forvalues j = 1/`n_vars_tr' {
			tempvar h_tmp_`j'
			qui gen double `h_tmp_`j'' = .
			local h_tmp_vars "`h_tmp_vars' `h_tmp_`j''"
		}

		tempname pipe_diag pipe_cvs pipe_tied
		matrix `pipe_diag' = J(4, 3, .)
		matrix `pipe_cvs' = J(4, 3, .)

		tempfile emat_file
		set seed `seed'
		_spur_backend write_emat, q(`q') nrep(`nrep') file(`"`emat_file'"')

		local latlong_flag = cond("`latlong'" != "", 1, 0)
		local plugin_vars "`slist' `varlist' `h_tmp_vars'"
		local pcmd `"spur_pipeline `d' `k_x_raw' 1 `latlong_flag' `q' `nrep' "`nthreads_arg'" `"`emat_file'"' "`pipe_diag'" "`pipe_cvs'" "`pipe_tied'""'
		capture noisily _spur_backend plugin_call, vars(`plugin_vars') cond(`touse_all') cmd(`"`pcmd'"')
		local rc = _rc
		capture erase `"`emat_file'"'
		if `rc' != 0 {
			exit `rc'
		}

		if "`req_backend'" == "auto" & scalar(`pipe_tied') == 1 {
			local use_pipeline 0
			local diag_tied 1
		}
		else {
			local diag_backend "rust"
			local diag_tied = scalar(`pipe_tied')

			scalar `i0_teststat' = `pipe_diag'[1, 1]
			scalar `i0_p' = `pipe_diag'[1, 2]
			scalar `i0_ha_param' = `pipe_diag'[1, 3]
			matrix `i0_cv' = (`pipe_cvs'[1, 1] \ `pipe_cvs'[1, 2] \ `pipe_cvs'[1, 3])

			scalar `i1_teststat' = `pipe_diag'[2, 1]
			scalar `i1_p' = `pipe_diag'[2, 2]
			scalar `i1_ha_param' = `pipe_diag'[2, 3]
			matrix `i1_cv' = (`pipe_cvs'[2, 1] \ `pipe_cvs'[2, 2] \ `pipe_cvs'[2, 3])

			scalar `i0resid_teststat' = `pipe_diag'[3, 1]
			scalar `i0resid_p' = `pipe_diag'[3, 2]
			scalar `i0resid_ha_param' = `pipe_diag'[3, 3]
			matrix `i0resid_cv' = (`pipe_cvs'[3, 1] \ `pipe_cvs'[3, 2] \ `pipe_cvs'[3, 3])

			scalar `i1resid_teststat' = `pipe_diag'[4, 1]
			scalar `i1resid_p' = `pipe_diag'[4, 2]
			scalar `i1resid_ha_param' = `pipe_diag'[4, 3]
			matrix `i1resid_cv' = (`pipe_cvs'[4, 1] \ `pipe_cvs'[4, 2] \ `pipe_cvs'[4, 3])

			matrix `diagnostics' = `pipe_diag'
			matrix rownames `diagnostics' = i0 i1 i0resid i1resid
			matrix colnames `diagnostics' = teststat p ha_param
		}
	}

	if `use_pipeline' == 0 {
		local indiv_backend = cond(`diag_tied' == 1, "fallback", "`req_backend'")
		local test_options q(`q') nrep(`nrep') seed(`seed') backend(`indiv_backend') nthreads(`nthreads_arg')
		if "`latlong'" != "" {
			local test_options `test_options' latlong
		}

		local transform_options prefix(`prefix') transformation(`transformation') backend(`indiv_backend') nthreads(`nthreads_arg')
		if "`latlong'" != "" {
			local transform_options `transform_options' latlong
		}
		if "`replace'" != "" {
			local transform_options `transform_options' replace
		}

		quietly spurtest i0 `depvar' `if' `in', `test_options'
		scalar `i0_teststat' = r(teststat)
		scalar `i0_p' = r(p)
		scalar `i0_ha_param' = r(ha_param)
		matrix `i0_cv' = r(cv)
		local be_i0 "`r(backend)'"
		local ti_i0 = r(tied)

		quietly spurtest i1 `depvar' `if' `in', `test_options'
		scalar `i1_teststat' = r(teststat)
		scalar `i1_p' = r(p)
		scalar `i1_ha_param' = r(ha_param)
		matrix `i1_cv' = r(cv)
		local be_i1 "`r(backend)'"
		local ti_i1 = r(tied)

		quietly spurtest i0resid `depvar' `indepvars' `if' `in', `test_options'
		scalar `i0resid_teststat' = r(teststat)
		scalar `i0resid_p' = r(p)
		scalar `i0resid_ha_param' = r(ha_param)
		matrix `i0resid_cv' = r(cv)
		local be_i0r "`r(backend)'"
		local ti_i0r = r(tied)

		quietly spurtest i1resid `depvar' `indepvars' `if' `in', `test_options'
		scalar `i1resid_teststat' = r(teststat)
		scalar `i1resid_p' = r(p)
		scalar `i1resid_ha_param' = r(ha_param)
		matrix `i1resid_cv' = r(cv)
		local be_i1r "`r(backend)'"
		local ti_i1r = r(tied)

		matrix `diagnostics' = ///
			(`i0_teststat', `i0_p', `i0_ha_param' \ ///
			 `i1_teststat', `i1_p', `i1_ha_param' \ ///
			 `i0resid_teststat', `i0resid_p', `i0resid_ha_param' \ ///
			 `i1resid_teststat', `i1resid_p', `i1resid_ha_param')
		matrix rownames `diagnostics' = i0 i1 i0resid i1resid
		matrix colnames `diagnostics' = teststat p ha_param

		if "`be_i0'" == "rust" & "`be_i1'" == "rust" & "`be_i0r'" == "rust" & "`be_i1r'" == "rust" {
			local diag_backend "rust"
		}
		else {
			local diag_backend "fallback"
		}
		local diag_tied = max(`diag_tied', `ti_i0', `ti_i1', `ti_i0r', `ti_i1r')
	}

	quietly regress `depvar' `indepvars' `if' `in', robust

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
	local levels_scpc_tied = e(tied)

	local transformed_depvar `prefix'`depvar'
	local transformed_indepvars
	foreach var of local indepvars {
		local transformed_indepvars `transformed_indepvars' `prefix'`var'
	}

	local tr_backend "`diag_backend'"
	if `use_pipeline' == 1 {
		forvalues j = 1/`n_vars_tr' {
			local vname : word `j' of `varlist'
			local hv : word `j' of `h_tmp_vars'
			if "`replace'" != "" {
				capture drop `prefix'`vname'
			}
			rename `hv' `prefix'`vname'
		}
	}
	else {
		quietly spurtransform `depvar' `indepvars' `if' `in', `transform_options'
		local tr_backend "`r(backend)'"
	}

	quietly regress `transformed_depvar' `transformed_indepvars' `if' `in', robust

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
	local tr_scpc_tied = e(tied)

	local used_backend "fallback"
	if "`diag_backend'" == "rust" & "`tr_backend'" == "rust" & "`levels_scpc_be'" == "rust" & "`tr_scpc_be'" == "rust" {
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

	return scalar i0_teststat = `i0_teststat'
	return scalar i0_p = `i0_p'
	return scalar i0_ha_param = `i0_ha_param'
	return matrix i0_cv = `i0_cv'

	return scalar i1_teststat = `i1_teststat'
	return scalar i1_p = `i1_p'
	return scalar i1_ha_param = `i1_ha_param'
	return matrix i1_cv = `i1_cv'

	return scalar i0resid_teststat = `i0resid_teststat'
	return scalar i0resid_p = `i0resid_p'
	return scalar i0resid_ha_param = `i0resid_ha_param'
	return matrix i0resid_cv = `i0resid_cv'

	return scalar i1resid_teststat = `i1resid_teststat'
	return scalar i1resid_p = `i1resid_p'
	return scalar i1resid_ha_param = `i1resid_ha_param'
	return matrix i1resid_cv = `i1resid_cv'

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
	return local scpc_levels_backend "`levels_scpc_be'"
	return local scpc_transformed_backend "`tr_scpc_be'"
	return local depvar "`depvar'"
	return local indepvars "`indepvars'"
	return local transformed_depvar "`transformed_depvar'"
	return local transformed_indepvars "`transformed_indepvars'"
	return local prefix "`prefix'"
	return local transformation "`transformation'"
end
