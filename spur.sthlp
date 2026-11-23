{smcl}
{* *! version 1.0.0 Apr 2026}{...}
{title:Title}

{pstd}
{hi:spur} {hline 2} Full SPUR pipeline wrapper.


{title:Syntax}

{p 8 16 2}
{cmd:spur} {depvar} {indepvars} {ifin} [{cmd:,} {it:options}]


{synoptset 13}{...}
{synopthdr}
{synoptline}
{synopt :{opt q(#)}}number of weighted averages used in the SPUR tests; default is {opt q(15)}{p_end}
{synopt :{opt nrep(#)}}number of Monte Carlo draws used in the SPUR tests; default is {opt nrep(100000)}{p_end}
{synopt :{opt seed(#)}}random seed used for the SPUR simulation draws; default is {opt seed(42)}{p_end}
{synopt :{opt latlong}}interpret {it:s_1} as latitude and {it:s_2} as longitude{p_end}
{synopt :{opt avc(#)}}maximal average pairwise correlation passed to {cmd:scpc}; default is {opt avc(0.03)}{p_end}
{synopt :{opt uncond}}pass {opt uncond} to {cmd:scpc}{p_end}
{synopt :{opt cvs}}pass {opt cvs} to {cmd:scpc} and store critical-value matrices{p_end}
{synopt :{opt replace}}allow overwriting transformed variables named {it:h_*}{p_end}
{synopt:{opt backend(auto|rust|fallback)}}choose the computational backend. {opt auto} uses the matching Rust plugin when available and Mata otherwise. {opt rust} requires the matching plugin. Set the default with {cmd:SPATIAL_BACKEND} as a Stata global or environment variable.{p_end}
{synopt:{opt nthreads(#)}}set the number of Rust threads. {opt nthreads(1)} runs Rust with one thread.{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:spur} runs the full SPUR workflow for a linear regression. It runs the
dependent-variable and residual SPUR diagnostics, estimates the regression in
levels, runs {cmd:scpc}, applies the default {cmd:lbmgls} transformation to the
dependent and independent variables, estimates the transformed regression, and
runs {cmd:scpc} again.

{pstd}
The wrapper uses observations with complete outcomes and regressors for all stages
and uses the same seeded simulation draws for all four diagnostics. It always
computes both the levels and transformed branches. Comparing both fits is a sensitivity analysis. Coverage after choosing a branch using the same data requires a separate justification.

{pstd}
The transformed variables are written to the active dataset with prefix
{it:h_}. Use {opt replace} if these variables already exist.
The transformed intercept reflects the centering of the LBM-GLS transformed variables.

{pstd}
{cmd:spur} requires the matching {cmd:scpc-stata} package. For IV and fixed-effect
workflows, fit the Stata model and call {cmd:scpc} directly.


{marker examples}{...}
{pstd}
SPUR diagnostics report the best simulated tuning power when the 50% target is out of reach. They identify nominal levels that the requested simulation count cannot attain. A fit using every positive direction or at least 40 components gives a high-frequency sensitivity advisory. Half-life inference flags limited information at 2–5 components and currently supports up to 100 components with its five-node quadrature.

{title:Examples}

{phang}{cmd:. spur am gini fracblack, latlong}{p_end}
{phang}{cmd:. spur am gini fracblack, q(10) nrep(10000) latlong replace}{p_end}
{phang}{cmd:. spur am gini fracblack, avc(0.01) uncond cvs latlong replace}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:spur} stores the following in {cmd:r()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 31 2: Scalars}{p_end}
{synopt:{cmd:r(tied)}}indicator of a retained tied eigenspace{p_end}
{synopt:{cmd:r(i0_teststat)}}I(0) diagnostic test statistic{p_end}
{synopt:{cmd:r(i0_p)}}I(0) diagnostic p-value{p_end}
{synopt:{cmd:r(i1_teststat)}}I(1) diagnostic test statistic{p_end}
{synopt:{cmd:r(i1_p)}}I(1) diagnostic p-value{p_end}
{synopt:{cmd:r(i0resid_teststat)}}residual I(0) diagnostic test statistic{p_end}
{synopt:{cmd:r(i0resid_p)}}residual I(0) diagnostic p-value{p_end}
{synopt:{cmd:r(i1resid_teststat)}}residual I(1) diagnostic test statistic{p_end}
{synopt:{cmd:r(i1resid_p)}}residual I(1) diagnostic p-value{p_end}
{synopt:{cmd:r(levels_N)}}number of observations in the levels regression{p_end}
{synopt:{cmd:r(levels_r2)}}R-squared from the levels regression{p_end}
{synopt:{cmd:r(transformed_N)}}number of observations in the transformed regression{p_end}
{synopt:{cmd:r(transformed_r2)}}R-squared from the transformed regression{p_end}

{p2col 5 28 31 2: Matrices}{p_end}
{synopt:{cmd:r(diagnostic_info)}}rows i0, i1, i0resid, i1resid and columns tied, all_directions, best_tuning_power{p_end}
{synopt:{cmd:r(diagnostics)}}diagnostic test statistics, p-values, and alternative parameters{p_end}
{synopt:{cmd:r(levels_scpcstats)}}SCPC inference table from the levels regression{p_end}
{synopt:{cmd:r(transformed_scpcstats)}}SCPC inference table from the transformed regression{p_end}
{synopt:{cmd:r(levels_scpccvs)}}levels critical values, only when {opt cvs} is specified{p_end}
{synopt:{cmd:r(transformed_scpccvs)}}transformed critical values, only when {opt cvs} is specified{p_end}

{p2col 5 28 31 2: Macros}{p_end}
{synopt:{cmd:r(backend)}}Rust if all stages used Rust, otherwise fallback{p_end}
{synopt:{cmd:r(diagnostics_backend)}}computational backend for the four diagnostics{p_end}
{synopt:{cmd:r(i0_backend)}}computational backend for the I(0) levels diagnostic{p_end}
{synopt:{cmd:r(i1_backend)}}computational backend for the I(1) levels diagnostic{p_end}
{synopt:{cmd:r(i0resid_backend)}}computational backend for the I(0) residual diagnostic{p_end}
{synopt:{cmd:r(i1resid_backend)}}computational backend for the I(1) residual diagnostic{p_end}
{synopt:{cmd:r(transform_backend)}}computational backend for the transformation{p_end}
{synopt:{cmd:r(scpc_levels_backend)}}computational backend for SCPC in levels{p_end}
{synopt:{cmd:r(scpc_transformed_backend)}}computational backend for transformed SCPC{p_end}
{synopt:{cmd:r(*_backend_reason)}}reason for using Mata for each stage, when applicable{p_end}
{synopt:{cmd:r(depvar)}}dependent variable{p_end}
{synopt:{cmd:r(indepvars)}}independent variables{p_end}
{synopt:{cmd:r(transformed_depvar)}}transformed dependent variable{p_end}
{synopt:{cmd:r(transformed_indepvars)}}transformed independent variables{p_end}
{synopt:{cmd:r(prefix)}}transformation prefix, currently {it:h_}{p_end}
{synopt:{cmd:r(transformation)}}transformation method, currently {it:lbmgls}{p_end}
{p2colreset}{...}


{marker authors}{...}
{title:Authors}

{pstd}
DGoettlich{break}
spatial-spur


{marker references}{...}
{title:References}

{phang}
Becker, Sascha O., P. David Boll and Hans-Joachim Voth "Testing and Correcting
for Spatial Unit Roots in Regression Analysis", The Stata Journal 26(2): 177–202. doi:10.1177/1536867X261449932. {browse "https://journals.sagepub.com/doi/10.1177/1536867X261449932"}.

{phang}
Müller, Ulrich K. and Mark W. Watson "Spatial Unit Roots and Spurious
Regression", Econometrica 92(5) (2024), 1661-1695. doi:10.3982/ECTA21654.
