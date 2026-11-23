*! version 1.1.0  17nov2026 spatial-spur
program define _spurtest_i0, rclass sortpreserve
	version 16.0
	syntax varname(numeric) [if] [in] , [q(int 15) nrep(int 100000) seed(int 42) latlong BACKend(string) nthreads(integer -1)]
	marksample touse

	_spur_backend run_test, flavor(i0) yvar(`varlist') touse(`touse') q(`q') nrep(`nrep') ///
		seed(`seed') `latlong' backend(`"`backend'"') nthreads(`nthreads')

	tempname p teststat g_a cv
	scalar `p' = r(p)
	scalar `teststat' = r(teststat)
	scalar `g_a' = r(ha_param)
	matrix `cv' = r(cv)
	local tied_val = r(tied)
	local all = r(all_directions)
	local best = r(best_tuning_power)
	local used_backend "`r(backend)'"
	local backend_reason "`r(backend_reason)'"

	return scalar p = `p'
	return scalar teststat = `teststat'
	return scalar ha_param = `g_a'
	return matrix cv = `cv'
	return scalar tied = `tied_val'
	return scalar all_directions = `all'
	return scalar best_tuning_power = `best'
	return local backend "`used_backend'"
	return local backend_reason "`backend_reason'"

	// Display results
	di as text "Spatial I(0) Test Results"
	di as result "---------------------------------------"
	di as result "Test Statistic (LFST) : " %9.4f `teststat'
	di as result "P-value               : " %9.4f `p'
	di as result "---------------------------------------"

end
