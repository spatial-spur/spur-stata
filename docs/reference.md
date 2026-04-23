# Reference

This page documents the public `spur-stata` command surface.

## Overview

| Command | Description |
|---|---|
| `spur` | Run the full SPUR workflow and return diagnostics plus levels and transformed branches |
| `spurtest` | Run one of the four SPUR diagnostic tests |
| `spurtransform` | Transform variables to remove spatial unit roots |
| `spurhalflife` | Estimate a spatial half-life confidence interval |

## Conventions

### Coordinates

All commands in the package expect the spatial coordinates to be stored in
variables named `s_*`.

- with `latlong`, `s_1` is latitude and `s_2` is longitude
- without `latlong`, the available `s_*` variables are treated as Euclidean
  coordinates

### Dependencies

- `moremata` is required by the SPUR commands
- `scpc` is a separate Stata package and handles the inference stage after the
  branch has been chosen

### Simulation Controls

The SPUR tests and half-life interval are simulation-based.

- `q(#)`: number of low-frequency weighted averages; default `q(15)`
- `nrep(#)`: number of Monte Carlo draws; default `nrep(100000)`

## Full Workflow

### `spur`

Run the full SPUR workflow for a linear regression.

**Syntax**

```stata
spur depvar indepvars [if] [in] [, q(#) nrep(#) latlong avc(#) uncond cvs replace]
```

**Behavior**

- runs `spurtest i0` and `spurtest i1` on the dependent variable
- runs `spurtest i0resid` and `spurtest i1resid` for the regression
- estimates the levels regression with `regress ..., robust`
- runs the external `scpc` postestimation command
- transforms the dependent and independent variables with `spurtransform`,
  using `prefix(h_)` and `transformation(lbmgls)`
- estimates the transformed regression with `regress ..., robust`
- runs `scpc` again

The wrapper always computes both the levels and transformed branches. It does
not encode the practitioner-guide branch decision as a returned field.

**Options**

- `q(#)`, `nrep(#)`, `latlong`: passed to the SPUR diagnostic tests
- `avc(#)`, `uncond`, `cvs`: passed to `scpc`
- `replace`: allow overwriting existing transformed `h_*` variables

**Stored results**

`spur` stores the following main results in `r()`:

- `r(diagnostics)`: rows `i0`, `i1`, `i0resid`, `i1resid`; columns
  `teststat`, `p`, `ha_param`
- `r(levels_scpcstats)`: SCPC inference table from the levels regression
- `r(transformed_scpcstats)`: SCPC inference table from the transformed
  regression
- `r(levels_scpccvs)` and `r(transformed_scpccvs)`: critical-value matrices
  when `cvs` is specified
- diagnostic scalar aliases such as `r(i0_p)` and `r(i1resid_teststat)`
- model-statistic scalars such as `r(levels_N)` and `r(transformed_r2)`

The wrapper is OLS-only. IV and fixed-effect workflows should be estimated
directly and followed by `scpc`.

## Diagnostics

### `spurtest`

Run one of the four SPUR tests.

**Syntax**

```stata
spurtest i1 varname [if] [in] [, q(#) nrep(#) latlong]
spurtest i0 varname [if] [in] [, q(#) nrep(#) latlong]
spurtest i1resid depvar [indepvars] [if] [in] [, q(#) nrep(#) latlong]
spurtest i0resid depvar [indepvars] [if] [in] [, q(#) nrep(#) latlong]
```

**Subcommands**

- `i0`: variable-level `I(0)` test
- `i1`: variable-level `I(1)` test
- `i0resid`: residual-based `I(0)` test for a linear regression
- `i1resid`: residual-based `I(1)` test for a linear regression

**Options**

- `q(#)`: number of weighted averages used in the test statistic
- `nrep(#)`: number of Monte Carlo draws used to simulate the reference
  distribution
- `latlong`: interpret `s_1` and `s_2` as latitude and longitude

**Stored results**

`spurtest` stores the following in `r()`:

- `r(teststat)`: LFUR or LFST test statistic
- `r(p)`: p-value
- `r(ha_param)`: alternative-hypothesis parameter (`c_a` or `g_a`)
- `r(cv)`: matrix of critical values at 1%, 5%, and 10%

The residual tests automatically include a constant in the underlying
regression.

## Transformation

### `spurtransform`

Transform one or more variables and write the transformed versions back into
the current dataset.

**Syntax**

```stata
spurtransform varlist [if] [in], prefix(string) ///
    [transformation(string) radius(#) clustvar(varname) latlong replace separately]
```

**Required option**

- `prefix(string)`: prefix applied to the transformed variables, for example
  `prefix(h_)` to create `h_am`

**Options**

- `transformation(string)`: transformation mode
  - `lbmgls`: GLS-style transformation; default branch
  - `nn`: nearest-neighbor differencing
  - `iso`: isotropic differencing within a radius
  - `cluster`: within-cluster demeaning
- `radius(#)`: radius for `transformation(iso)`; required there and not allowed
  otherwise
- `clustvar(varname)`: cluster variable for `transformation(cluster)`; required
  there and not allowed otherwise
- `latlong`: interpret `s_1` and `s_2` as latitude and longitude
- `replace`: allow overwriting previously created transformed variables
- `separately`: construct the transformation separately for each variable
  rather than on the common non-missing sample

**Behavior**

- the command transforms every variable in `varlist`
- transformed variables are written back into the active dataset
- there are no stored results; the output is the newly created variables

## Persistence

### `spurhalflife`

Estimate a confidence interval for the spatial half-life of a variable.

**Syntax**

```stata
spurhalflife varname [if] [in] [, q(#) nrep(#) level(#) latlong normdist]
```

**Options**

- `q(#)`: number of weighted averages; default `q(15)`
- `nrep(#)`: number of Monte Carlo draws; default `nrep(100000)`
- `level(#)`: confidence level in percent; default `level(95)`
- `latlong`: interpret `s_1` and `s_2` as latitude and longitude
- `normdist`: report the interval as a fraction of the maximum pairwise
  distance rather than in distance units

**Stored results**

`spurhalflife` stores the following in `r()`:

- `r(ci_l)`: lower interval bound
- `r(ci_u)`: upper interval bound
- `r(max_dist)`: maximum pairwise distance in the sample

If `normdist` is not specified, the interval is reported in metres under
`latlong`, or in the units of the Euclidean coordinates otherwise.

## Workflow Note

Use `spur` for the full OLS pipeline in one command. Use `spurtest`,
`spurtransform`, `spurhalflife`, and direct `scpc` calls when you want only
selected steps or when you need IV / fixed-effect estimation.
