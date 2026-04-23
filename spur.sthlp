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
{synopt :{opt latlong}}interpret {it:s_1} as latitude and {it:s_2} as longitude{p_end}
{synopt :{opt avc(#)}}maximal average pairwise correlation passed to {cmd:scpc}; default is {opt avc(0.03)}{p_end}
{synopt :{opt uncond}}pass {opt uncond} to {cmd:scpc}{p_end}
{synopt :{opt cvs}}pass {opt cvs} to {cmd:scpc} and store critical-value matrices{p_end}
{synopt :{opt replace}}allow overwriting transformed variables named {it:h_*}{p_end}
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
The wrapper is intended to mirror the R and Python {cmd:spur()} wrappers. It
always computes both the levels and transformed branches. The practitioner
decision rule is not encoded as a returned field.

{pstd}
The transformed variables are written to the active dataset with prefix
{it:h_}. Use {opt replace} if these variables already exist.

{pstd}
{cmd:scpc} is an external postestimation command from the Müller-Watson
{cmd:SCPC} package. IV and fixed-effect workflows should be estimated directly
and followed by {cmd:scpc}; they are not handled by this wrapper.


{marker examples}{...}
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
{synopt:{cmd:r(diagnostics)}}diagnostic test statistics, p-values, and alternative parameters{p_end}
{synopt:{cmd:r(levels_scpcstats)}}SCPC inference table from the levels regression{p_end}
{synopt:{cmd:r(transformed_scpcstats)}}SCPC inference table from the transformed regression{p_end}
{synopt:{cmd:r(levels_scpccvs)}}levels critical values, only when {opt cvs} is specified{p_end}
{synopt:{cmd:r(transformed_scpccvs)}}transformed critical values, only when {opt cvs} is specified{p_end}

{p2col 5 28 31 2: Macros}{p_end}
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
for Spatial Unit Roots in Regression Analysis", Stata Journal, forthcoming.

{phang}
Müller, Ulrich K. and Mark W. Watson "Spatial Unit Roots and Spurious
Regression", Econometrica 92 (2024), 1661-1695.
