*! version 1.1.0  17nov2026 spatial-spur
mata mata set matastrict on

program define spurhalflife, rclass sortpreserve
	version 14
	syntax varname(numeric) [if] [in] , [q(int 15) nrep(int 100000) seed(int 42) Level(real 95) latlong NORMdist BACKend(string) nthreads(integer -1)]
	marksample touse

	if `level' >= 100 | `level' <= 0 {
		di as error "Invalid level."
		exit 198
	}
	if `q' < 1 {
		di as error "q() must be at least 1"
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

	_spur_backend resolve, backend(`"`backend'"') nthreads(`nthreads')
	local req_backend "`r(req_backend)'"
	local nthreads_arg "`r(nthreads)'"
	local plugin_ok = r(plugin_ok)

	_spur_backend check_coords, touse(`touse') `latlong' quiet
	local slist "`r(slist)'"
	local d = r(d)
	local n_obs = r(n)

	tempname ci_l ci_u max_dist tied_scal
	local used_backend ""
	local tied_val 0

	if "`req_backend'" == "fallback" | `plugin_ok' == 0 {
		_spur_backend check_moremata
		set seed `seed'
		_spur_backend mata_halflife, yvar(`varlist') touse(`touse') `latlong' q(`q') nrep(`nrep') ///
			level(`level') cil(`ci_l') ciu(`ci_u') maxdist(`max_dist')
		local used_backend "fallback"
	}
	else {
		tempfile emat_file
		set seed `seed'
		_spur_backend write_emat, q(`q') nrep(`nrep') file(`"`emat_file'"')

		local latlong_flag = cond("`latlong'" != "", 1, 0)
		local pcmd `"spur_halflife `d' `latlong_flag' `q' `nrep' "`level'" "`nthreads_arg'" `"`emat_file'"' "`ci_l'" "`ci_u'" "`max_dist'" "`tied_scal'""'
		capture noisily _spur_backend plugin_call, vars(`slist' `varlist') cond(`touse') cmd(`"`pcmd'"')
		local rc = _rc
		capture erase `"`emat_file'"'
		if `rc' != 0 {
			exit `rc'
		}

		if "`req_backend'" == "auto" & scalar(`tied_scal') == 1 {
			_spur_backend check_moremata
			set seed `seed'
			_spur_backend mata_halflife, yvar(`varlist') touse(`touse') `latlong' q(`q') nrep(`nrep') ///
				level(`level') cil(`ci_l') ciu(`ci_u') maxdist(`max_dist')
			local used_backend "fallback"
			local tied_val 1
		}
		else {
			_spur_backend print_coords, n(`n_obs') d(`d') `latlong'
			local used_backend "rust"
			local tied_val = scalar(`tied_scal')
		}
	}

	if !missing(scalar(`ci_u')) & scalar(`ci_u') >= 100 {
		scalar `ci_u' = .
	}

	if "`normdist'" == "" {
		scalar `ci_l' = scalar(`ci_l') * scalar(`max_dist')
		scalar `ci_u' = scalar(`ci_u') * scalar(`max_dist')
	}

	return scalar ci_l = scalar(`ci_l')
	return scalar ci_u = scalar(`ci_u')
	return scalar max_dist = scalar(`max_dist')
	return scalar tied = `tied_val'
	return local backend "`used_backend'"

	// Display results
	if "`normdist'" != "" {
		local units : di "fractions of maximum distance " %9.4f scalar(`max_dist')
	}
	else {
		if "`latlong'" != "" {
			local units "metres"
		}
		else {
			local units "unit of original coordinates"
		}
	}
	di as text "Spatial half life " "`level'" "% confidence interval (in " "`units'" ")"
	di as result "---------------------------------------"
	di as result "Lower bound: " %9.4f scalar(`ci_l')
	if !missing(scalar(`ci_u')) {
		di as result "Upper bound: " %9.4f scalar(`ci_u')
	}
	else {
		di as result "Upper bound:    inf"
	}
	di as result "---------------------------------------"

end
