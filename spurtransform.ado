*! version 1.1.0  17nov2026 spatial-spur
mata mata set matastrict on

program spurtransform, rclass sortpreserve
	version 14
	syntax varlist(numeric) [if] [in] , PREfix(string) [ transformation(string) radius(real -1) clustvar(varname) latlong Replace separately BACKend(string) nthreads(integer -1)]

	// check correct specification of options
	if "`transformation'" == "" {
		local transformation "lbmgls"
	}

	if "`transformation'" != "lbmgls" & "`transformation'" != "nn" & "`transformation'" != "iso" & "`transformation'" != "cluster" {
		di as error "Invalid transformation: `transformation' (expected lbmgls, nn, iso, or cluster)"
		exit 198
	}

	if "`transformation'" != "iso" & (`radius' != -1 | missing(`radius')) {
		di as error "Option radius only allowed with transformation(iso)."
		exit 198
	}

	if "`transformation'" == "iso" & `radius' == -1 {
		di as error "Radius missing."
		exit 198
	}

	if "`transformation'" == "iso" & (`radius' <= 0 | missing(`radius')) {
		di as error "Radius must be positive."
		exit 198
	}

	if "`transformation'" != "cluster" & "`clustvar'" != "" {
		di as error "Option clustvar only allowed with transformation(cluster)."
		exit 198
	}

	if "`transformation'" == "cluster" & "`clustvar'" == "" {
		di as error "Clustvar missing."
		exit 198
	}

	_spur_backend resolve, backend(`"`backend'"') nthreads(`nthreads')
	local req_backend "`r(req_backend)'"
	local nthreads_arg "`r(nthreads)'"
	local plugin_ok = r(plugin_ok)
	local backend_reason "`r(backend_reason)'"

	// make numeric clustervar
	tempvar clustervar
	local has_cluster 0
	if "`transformation'" == "cluster" & "`clustvar'" != "" {
		qui egen `clustervar' = group(`clustvar')
		local has_cluster 1
	}

	// select sample
	if "`separately'" != "" {
		marksample touse, novarlist
	}
	else {
		marksample touse
	}

	confirm name `prefix'
	foreach name of local varlist {
		capture confirm new variable `prefix'`name'
		if _rc {
			if "`replace'" == "" {
				di as error "`prefix'`name' already defined or invalid name"
				exit 110
			}
		}
	}

	local used_backend ""
	if "`req_backend'" == "fallback" | `plugin_ok' == 0 {
		local used_backend "fallback"
	}
	else {
		local used_backend "rust"
	}

	local latlong_flag = cond("`latlong'" != "", 1, 0)
	local cl_arg = cond(`has_cluster' == 1, "`clustervar'", "")

	local n_y : word count `varlist'
	local out_vars ""
	forvalues j = 1/`n_y' {
		tempvar out_`j'
		qui gen double `out_`j'' = .
		local out_vars "`out_vars' `out_`j''"
	}

	tempvar touse2
	if "`separately'" == "" {
		qui gen byte `touse2' = `touse'
		if "`transformation'" == "cluster" {
			qui replace `touse2' = 0 if missing(`clustervar')
		}

		_spur_backend check_coords, touse(`touse2') `latlong' quiet
		local slist "`r(slist)'"
		local d = r(d)
		local n_obs = r(n)

		if "`used_backend'" == "rust" {
			_spur_backend print_coords, n(`n_obs') d(`d') `latlong'
			local pcmd `"spur_transform `d' "`transformation'" `latlong_flag' "`radius'" `has_cluster' `n_y' "`nthreads_arg'""'
			_spur_backend plugin_call, vars(`slist' `cl_arg' `varlist' `out_vars') cond(`touse2') cmd(`"`pcmd'"')
		}
		else {
			_spur_backend mata_transform, invars(`varlist') outvars(`out_vars') touse(`touse2') ///
				transformation(`transformation') radius(`radius') clustvar(`cl_arg') `latlong'
		}
	}
	else {
		forvalues j = 1/`n_y' {
			local name : word `j' of `varlist'
			local ov : word `j' of `out_vars'
			qui gen byte `touse2' = `touse'
			qui replace `touse2' = 0 if missing(`name')
			if "`transformation'" == "cluster" {
				qui replace `touse2' = 0 if missing(`clustervar')
			}

			_spur_backend check_coords, touse(`touse2') `latlong' quiet
			local slist "`r(slist)'"
			local d = r(d)
			local n_obs = r(n)

			if "`used_backend'" == "rust" {
				_spur_backend print_coords, n(`n_obs') d(`d') `latlong'
				local pcmd `"spur_transform `d' "`transformation'" `latlong_flag' "`radius'" `has_cluster' 1 "`nthreads_arg'""'
				_spur_backend plugin_call, vars(`slist' `cl_arg' `name' `ov') cond(`touse2') cmd(`"`pcmd'"')
			}
			else {
				_spur_backend mata_transform, invars(`name') outvars(`ov') touse(`touse2') ///
					transformation(`transformation') radius(`radius') clustvar(`cl_arg') `latlong'
			}
			qui drop `touse2'
		}
	}

	// Commit outputs only after all computations succeed
	forvalues j = 1/`n_y' {
		local name : word `j' of `varlist'
		local ov : word `j' of `out_vars'
		if "`replace'" != "" {
			capture drop `prefix'`name'
		}
		rename `ov' `prefix'`name'
	}

	return local backend "`used_backend'"
	return local backend_reason "`backend_reason'"
	return local transformation "`transformation'"
end
