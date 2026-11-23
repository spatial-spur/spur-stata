*! version 1.1.0  17nov2026 spatial-spur
mata mata set matastrict on

program define spurhalflife, rclass sortpreserve
	version 16.0
	syntax varname(numeric) [if] [in] , [q(int 15) nrep(int 100000) seed(int 42) Level(real 95) latlong NORMdist BACKend(string) nthreads(integer -1)]
	marksample touse

	if `level' >= 100 | `level' <= 0 {
		di as error "Invalid level."
		exit 198
	}
	if `q' < 2 {
		di as error "q=`q' must be at least 2 for half-life inference"
		exit 198
	}
    if `q' > 100 {
        di as error "q=`q' exceeds the validated half-life quadrature limit of 100"
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
	local backend_reason "`r(backend_reason)'"

	_spur_backend check_coords, touse(`touse') `latlong' quiet minobs(3)
	local slist "`r(slist)'"
	local d = r(d)
	local n_obs = r(n)

	if `q' >= `n_obs' {
		di as error "q=`q' must be smaller than the number of observations (n=`n_obs')"
		exit 198
	}

	tempname ci_l ci_u max_dist info components ncomp
	local used_backend ""
	matrix `info' = J(1,3,.)

	if "`req_backend'" == "fallback" | `plugin_ok' == 0 {
		set seed `seed'
		_spur_backend mata_halflife, yvar(`varlist') touse(`touse') `latlong' q(`q') nrep(`nrep') ///
			level(`level') cil(`ci_l') ciu(`ci_u') maxdist(`max_dist') components(`components') ncomp(`ncomp') infomat(`info')
		local used_backend "fallback"
	}
	else {
		matrix `components' = J(65, 2, .)
		tempfile emat_file
		set seed `seed'
		_spur_backend write_emat, q(`q') nrep(`nrep') file(`"`emat_file'"')

		local latlong_flag = cond("`latlong'" != "", 1, 0)
		local pcmd `"spur_halflife `d' `latlong_flag' `q' `nrep' "`level'" "`nthreads_arg'" `"`emat_file'"' "`ci_l'" "`ci_u'" "`max_dist'" "`info'" "`components'" "`ncomp'""'
		capture noisily _spur_backend plugin_call, vars(`slist' `varlist') cond(`touse') cmd(`"`pcmd'"')
		local rc = _rc
		capture erase `"`emat_file'"'
		if `rc' != 0 {
			exit `rc'
		}

        _spur_backend print_coords, n(`n_obs') d(`d') `latlong'
        local used_backend "rust"
    }

    local n_components = scalar(`ncomp')
    local empty = `n_components' == 0
    local unbounded = !`empty' & missing(scalar(`ci_u'))
    if !`empty' {
        matrix `components' = `components'[1..`n_components', 1..2]
        if "`normdist'" == "" {
            mata: st_matrix("`components'", st_matrix("`components'") * st_numscalar("`max_dist'"))
        }
        matrix colnames `components' = lower upper
    }

	if "`normdist'" == "" {
		scalar `ci_l' = scalar(`ci_l') * scalar(`max_dist')
		scalar `ci_u' = scalar(`ci_u') * scalar(`max_dist')
	}

	_spur_backend inference_note, q(`q') nrep(`nrep') info(`info') level(`level')
	return scalar ci_l = scalar(`ci_l')
	return scalar ci_u = scalar(`ci_u')
	return scalar max_dist = scalar(`max_dist')
	return scalar empty = `empty'
    return scalar unbounded = `unbounded'
    return scalar n_components = `n_components'
	return scalar tied = `info'[1,1]
	return scalar all_directions = `info'[1,2]
	return local backend "`used_backend'"
	return local backend_reason "`backend_reason'"

	// Display results
	if "`normdist'" != "" {
		local max_str = strtrim("`: di %12.6g scalar(`max_dist')'")
		local units "fractions of maximum distance `max_str'"
	}
	else {
		if "`latlong'" != "" {
			local units "metres"
		}
		else {
			local units "unit of original coordinates"
		}
	}
    di as text "Approximate spatial half-life `level'% confidence set (in `units')"
    di as result "---------------------------------------"
    if `empty' {
        di as result "Empty confidence set."
    }
    else {
        forvalues i = 1/`n_components' {
            local upper = cond(missing(`components'[`i', 2]), "Inf", strofreal(`components'[`i', 2], "%12.6g"))
            di as result "[" %12.6g `components'[`i', 1] ", " strtrim("`upper'") "]"
        }
    }
    di as result "---------------------------------------"
    if !`empty' {
        return matrix components = `components'
    }

end
