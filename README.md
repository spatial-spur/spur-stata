<p align="center">
  <img src="assets/logo.png" alt="SPUR logo">
</p>

# `spur`: A Stata Package around Spatial Unit Roots

This package implements methods for diagnosing and correcting spatial unit roots developed by Müller & Watson 2024. A practical guide to these methods and the Stata implementation can be found in [Becker, Boll and Voth 2026](https://journals.sagepub.com/doi/10.1177/1536867X261449932).

**Citation:** If you use this package, please cite Becker, Boll and Voth 2026 and Müller & Watson 2024. See [CITATION.bib](CITATION.bib) for the BibTeX entries.

If you encounter any issues or have any questions, please open an issue on GitHub or contact the authors.

## Installation

Installation from GitHub:

    . net install spur, replace from(https://raw.githubusercontent.com/spatial-spur/spur-stata/v0.1.1/)

If you *also* want to download the example data and do-file (will be saved to the current directory), run the following *in addition to the above*:

    . net get spur, replace from(https://raw.githubusercontent.com/spatial-spur/spur-stata/v0.1.1/)

Then, you can run the example do-file, which reproduces Table 1 from [Müller & Watson 2024](https://www.princeton.edu/~umueller/SPUR.pdf) based on the data from [Chetty et. al. (2014)](https://doi.org/10.1093/qje/qju022) (see [Becker, Boll and Voth 2026](https://journals.sagepub.com/doi/10.1177/1536867X261449932) for details):

    . do example

### Local development checkout

Install the matching SCPC checkout before using the combined `spur` command.

```stata
net install scpc, replace from("/path/to/spatial-spur/scpc-stata/src")
net install spur, replace from("/path/to/spatial-spur/spur-stata")
```

The Rust plugins use OxiBLAS through the shared Rust core.
The [OxiBLAS license](spur-oxiblas-license.txt) is included in the installation.

The development checkout bundles plugins for Apple Silicon, tested with
StataNow 19.5. `net install` installs the plugins on Apple Silicon and the
Mata implementations on all platforms. Build both plugins with
`make -C spatial-rust/crates/spatial-stata` from the shared workspace
after fetching the locked Cargo dependencies. The Makefile has macOS and Linux
targets. OxiBLAS selects SME2 at runtime when its hardware and operating-system
checks pass. Set `OXIBLAS_DISABLE_SME2=1` in the environment before starting Stata
to opt out. Use `adopath ++` with the `scpc-stata/src` and `spur-stata` directories
to load locally built plugins from your checkouts, then run
`mata: mata mlib index` to index the SPUR Mata library. To rebuild that library,
run `stata-se -b do mata/build_lspur.do` from the SPUR checkout in a fresh session.

Restart Stata after replacing a plugin. The default `backend(auto)` uses Mata
if a matching plugin cannot load. You can select Mata with `backend(fallback)`.
Check `r(backend)` and `r(backend_reason)` after standalone commands. The
combined `spur` command records the backend and fallback reason for each stage.

## Basic Usage

The package includes one wrapper `spur` implementing the whole workflow, and the individual commands `spurtest`, `spurtransform`, and `spurhalflife`. Example uses of these commands:.

```stata
// Run the full OLS workflow
spur am gini fracblack, latlong replace

// Test the null hypothesis of a spatial unit root in the variable am
spurtest i1 am, latlong

// Transform the variable am to remove the spatial unit root and save the transformed variable as h_am
spurtransform am, prefix("h_") latlong

// Construct a confindence interval for the spatial half-life of the transformed variable am
spurhalflife am, latlong                

```
For details, see the help files, as well as [Becker, Boll and Voth 2026](https://journals.sagepub.com/doi/10.1177/1536867X261449932).

## References

Becker, Sascha O., P. David Boll and Hans-Joachim Voth "Testing and Correcting for Spatial Unit Roots in Regression Analysis", The Stata Journal 26(2): 177–202. [https://doi.org/10.1177/1536867X261449932](https://journals.sagepub.com/doi/10.1177/1536867X261449932).

Chetty, Raj, Nathaniel Hendren, Patrick Kline, Emmanuel Saez "Where is the land of Opportunity? The Geography of Intergenerational Mobility in the United States" , The Quarterly Journal of Economics 129(4) (2014), 1553–1623, https://doi.org/10.1093/qje/qju022

Müller, Ulrich K. and Mark W. Watson "Spatial Unit Roots and Spurious Regression", Econometrica 92(5) (2024), 1661–1695. https://www.princeton.edu/~umueller/SPUR.pdf.

SPUR inference requires at least two informative weighted directions. Residual
procedures also require `q` no larger than the residual degrees of freedom.
Experiments with insufficient covariance rank or proportional covariance shapes
return an error. Euclidean half-lives use the coordinate units; geographic
half-lives use meters. The Mata fallback uses the same type-7 simulation
quantiles and canonical basis orientation as the native backend; these choices
can change finite-simulation results relative to earlier Mata releases.
