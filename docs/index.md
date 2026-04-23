# spur-stata

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

This walkthrough follows the practitioner guide in Becker, Boll, and Voth
(2026). The branch decision is based on the dependent-variable `I(0)` and
`I(1)` tests, using a 10% significance level.

## 1. Prepare the sample

Construct the estimation sample and retain the variables used in the
diagnostics and regression.

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

`spur-stata` expects the coordinate variables to be named `s_1`, `s_2`. 
With `latlong`, `s_1` is latitude and `s_2` is longitude.

## 2. Run the full pipeline wrapper

Use `spur` when you want the full OLS workflow in one command. It runs the
diagnostics, estimates the levels branch, runs `scpc`, transforms the model
variables with `lbmgls`, estimates the transformed branch, and runs `scpc`
again.

```stata
spur am gini fracblack, q(15) nrep(100000) latlong replace
```

The wrapper stores diagnostics in `r(diagnostics)`, levels inference in
`r(levels_scpcstats)`, and transformed inference in
`r(transformed_scpcstats)`.

## 3. Run individual commands instead

Use the individual commands when you want only selected parts of the workflow.

### Test the dependent variable against the I(0) alternative

The first diagnostic tests the dependent variable under the spatial `I(0)`
null.

```stata
spurtest i0 am, latlong
```

The command stores the main results in `r()`, including `r(teststat)` and
`r(p)`.

### Test the dependent variable against the I(1) alternative

The second diagnostic tests the same variable under the spatial `I(1)` null.

```stata
spurtest i1 am, latlong
```

### Apply the decision rule

Using a 10% significance threshold:

- If you do **not** reject `I(0)` and you **do** reject `I(1)`, proceed in
  levels.
- In every other case, treat the specification as requiring spatial
  differencing and transform the dependent and independent variables together.

### Levels branch

Use this branch only when the decision rule implies that the dependent variable
is consistent with `I(0)` and inconsistent with `I(1)`.

```stata
reg am gini fracblack, robust
scpc, latlong
```

In this branch there is no SPUR transformation step; the regression is
estimated in levels and `scpc` is used for inference.

### Transformed branch

In every other case, transform the dependent and independent variables
together, re-estimate the regression on the transformed data, and use `scpc`
for inference there.

```stata
spurtransform am gini fracblack, prefix(h_) transformation(lbmgls) latlong replace

reg h_am h_gini h_fracblack, robust
scpc, latlong
```

The default empirical branch is the `lbmgls` transformation.

### Residual diagnostics

`spur-stata` also provides residual-based `I(0)` and `I(1)` tests:

```stata
spurtest i0resid am gini fracblack, latlong
spurtest i1resid am gini fracblack, latlong
```

## 4. SCPC-only workflows

The `spur` wrapper is OLS-only. For IV or fixed-effect models, estimate the
model directly with the relevant Stata command and run `scpc` immediately
afterward.

## Next Step

See [Reference](reference.md) for the full command syntax, option meanings, and
stored results.
