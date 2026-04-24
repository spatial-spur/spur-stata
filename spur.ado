*! version 1.0.0  23apr2026 DGoettlich
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
		]

	// check moremata installed
	cap which moremata.hlp
	if _rc {
		display as error in smcl `"Please install package {it:moremata} from SSC in order to use this command;"' _newline ///
        `"you can do so by clicking this link: {stata "ssc install moremata, replace":auto-install moremata}"'
		exit 199
	}

	// check scpc installed
	cap which scpc
	if _rc {
		display as error in smcl `"Please install package {it:scpc} in order to use this command;"' _newline ///
        `"which you can find under: {browse "https://github.com/ukmueller/SCPC/"}"'
		exit 199
	}

	if `avc'==-1 {
		local avc 0.03
	}

	gettoken depvar indepvars : varlist

	local prefix h_
	local transformation lbmgls
	local ncoef : word count `varlist'

	local test_options q(`q') nrep(`nrep') seed(`seed')
	if "`latlong'"!="" {
		local test_options `test_options' latlong
	}

	local scpc_options avc(`avc') k(`ncoef')
	if "`latlong'"!="" {
		local scpc_options `scpc_options' latlong
	}
	if "`uncond'"!="" {
		local scpc_options `scpc_options' uncond
	}
	if "`cvs'"!="" {
		local scpc_options `scpc_options' cvs
	}

	local transform_options prefix(`prefix') transformation(`transformation')
	if "`latlong'"!="" {
		local transform_options `transform_options' latlong
	}
	if "`replace'"!="" {
		local transform_options `transform_options' replace
	}

	tempname i0_teststat i0_p i0_ha_param i0_cv
	quietly spurtest i0 `depvar' `if' `in', `test_options'
	scalar `i0_teststat' = r(teststat)
	scalar `i0_p' = r(p)
	scalar `i0_ha_param' = r(ha_param)
	matrix `i0_cv' = r(cv)

	tempname i1_teststat i1_p i1_ha_param i1_cv
	quietly spurtest i1 `depvar' `if' `in', `test_options'
	scalar `i1_teststat' = r(teststat)
	scalar `i1_p' = r(p)
	scalar `i1_ha_param' = r(ha_param)
	matrix `i1_cv' = r(cv)

	tempname i0resid_teststat i0resid_p i0resid_ha_param i0resid_cv
	quietly spurtest i0resid `depvar' `indepvars' `if' `in', `test_options'
	scalar `i0resid_teststat' = r(teststat)
	scalar `i0resid_p' = r(p)
	scalar `i0resid_ha_param' = r(ha_param)
	matrix `i0resid_cv' = r(cv)

	tempname i1resid_teststat i1resid_p i1resid_ha_param i1resid_cv
	quietly spurtest i1resid `depvar' `indepvars' `if' `in', `test_options'
	scalar `i1resid_teststat' = r(teststat)
	scalar `i1resid_p' = r(p)
	scalar `i1resid_ha_param' = r(ha_param)
	matrix `i1resid_cv' = r(cv)

	tempname diagnostics
	matrix `diagnostics' = ///
		(`i0_teststat', `i0_p', `i0_ha_param' \ ///
		 `i1_teststat', `i1_p', `i1_ha_param' \ ///
		 `i0resid_teststat', `i0resid_p', `i0resid_ha_param' \ ///
		 `i1resid_teststat', `i1resid_p', `i1resid_ha_param')
	matrix rownames `diagnostics' = i0 i1 i0resid i1resid
	matrix colnames `diagnostics' = teststat p ha_param

	quietly regress `depvar' `indepvars' `if' `in', robust

	tempname levels_N levels_r2 levels_r2_a
	scalar `levels_N' = e(N)
	scalar `levels_r2' = e(r2)
	scalar `levels_r2_a' = e(r2_a)

	quietly scpc, `scpc_options'

	tempname levels_scpcstats levels_scpccvs
	matrix `levels_scpcstats' = e(scpcstats)
	if "`cvs'"!="" {
		matrix `levels_scpccvs' = e(scpccvs)
	}

	quietly spurtransform `depvar' `indepvars' `if' `in', `transform_options'

	local transformed_depvar `prefix'`depvar'
	local transformed_indepvars
	foreach var of local indepvars {
		local transformed_indepvars `transformed_indepvars' `prefix'`var'
	}

	quietly regress `transformed_depvar' `transformed_indepvars' `if' `in', robust

	tempname transformed_N transformed_r2 transformed_r2_a
	scalar `transformed_N' = e(N)
	scalar `transformed_r2' = e(r2)
	scalar `transformed_r2_a' = e(r2_a)

	quietly scpc, `scpc_options'

	tempname transformed_scpcstats transformed_scpccvs
	matrix `transformed_scpcstats' = e(scpcstats)
	if "`cvs'"!="" {
		matrix `transformed_scpccvs' = e(scpccvs)
	}

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

	if "`cvs'"!="" {
		return matrix levels_scpccvs = `levels_scpccvs'
		return matrix transformed_scpccvs = `transformed_scpccvs'
	}

	return local depvar "`depvar'"
	return local indepvars "`indepvars'"
	return local transformed_depvar "`transformed_depvar'"
	return local transformed_indepvars "`transformed_indepvars'"
	return local prefix "`prefix'"
	return local transformation "`transformation'"
end
