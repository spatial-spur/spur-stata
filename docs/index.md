<div class="package-overview-header" markdown="1">

# spur-stata

[GitHub repository](https://github.com/spatial-spur/spur-stata){ .md-button }

</div>

`spur-stata` implements the SPUR workflow for diagnosing and correcting
spatial unit roots in cross-sectional regressions. It provides both a full OLS
pipeline wrapper and the individual diagnostic, transformation, and half-life
commands. For inference, it is designed to work alongside the `scpc` package by
Müller and Watson.

## Installation

```stata
ssc install moremata, replace
net install spur, replace from(https://raw.githubusercontent.com/spatial-spur/spur-stata/v0.1.0b1/)

cap ado uninstall scpc
net install scpc, from("https://raw.githubusercontent.com/ukmueller/SCPC/master/src")
```

## Example: Chetty Dataset

In this example, we walk you through the workflow we recommend with the
packages step-by-step. We also provide a one-stop [pipeline wrapper](#pipeline-wrapper)
implementing the entire workflow in one step.

### Data preparation

For illustration, we load the Chetty dataset we ship as part of the package. Of
course, the analysis in principle follows the same logic on any other dataset.
In this specific case, we first omit the non-contiguous US states. We also drop
rows with missing values.

```stata
import excel "example_data/Chetty_Data_1.xlsx", ///
    sheet("Sheet1") firstrow case(lower) clear

drop if state == "HI"
drop if state == "AK"

rename lat s_1
rename lon s_2

keep am gini fracblack s_1 s_2 state
drop if missing(am, gini, fracblack, s_1, s_2)
```

In Stata, the SPUR and SCPC commands expect coordinate variables named `s_1`
and `s_2`.
With `latlong`, `s_1` is latitude and `s_2` is longitude.

### Testing for a spatial unit root

Based on MW 2024, we suggest first testing for a spatial unit root setting
using the `I(0)` and `I(1)` tests on the dependent variable.

One way to do this is to use the `spurtest i0` and `spurtest i1` commands
directly:

```stata
// am is the dependent variable
spurtest i0 am, latlong
spurtest i1 am, latlong
```

### Interpreting the test statistics

Using a 10% significance threshold, we suggest interpreting the results with the following heuristic:

- If you do **not** reject `I(0)` and you **do** reject `I(1)`, there is **likely no spatial unit root** and you can proceed in levels
- every other case implies a **possible spatial unit root** - in that case, we suggest transforming all dependent and independent variables before running regressions

We suggest always applying SCPC inference.

### Case 1: likely no spatial unit root

If the heuristic implies your scenario is unlikely to be a spatial unit root, we suggest proceeding in levels but applying SCPC inference:

```stata
regress am gini fracblack, robust
scpc, latlong
```

### Case 2: likely spatial unit root

If you do have a likely spatial unit root according to the heuristic above, we suggest applying the transformation and running the regression on transformed variables with SCPC inference:

```stata
spurtransform am gini fracblack, prefix(h_) transformation(lbmgls) latlong replace

regress h_am h_gini h_fracblack, robust
scpc, latlong
```

### Pipeline wrapper

As a shortcut to implementing all of those steps individually, we also provide a
`spur` wrapper that implements the entire pipeline in one step. It simply runs
all tests and returns all results.

```stata
spur am gini fracblack, latlong replace
```

The wrapper stores diagnostics in `r(diagnostics)`, levels inference in
`r(levels_scpcstats)`, and transformed inference in
`r(transformed_scpcstats)`.

### Residual tests

We also provide tests for spatial unit roots in regression residuals
rather than the dependent variable itself:

```stata
spurtest i0resid am gini fracblack, latlong
spurtest i1resid am gini fracblack, latlong
```

## SCPC-only workflows

The `spur` wrapper is OLS-only. For IV models, estimate the
model directly with the relevant Stata command and run `scpc` immediately
afterward.

## Next Step

See [Reference](reference.md) for the full command syntax, option meanings, and
stored results.
