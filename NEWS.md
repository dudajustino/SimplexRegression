# SimplexRegression 0.1.2

* Removed redundant "Provides functions for" from DESCRIPTION.
* Added references to DESCRIPTION (Barndorff-Nielsen & Jorgensen, 1991;
  Justino & Cribari-Neto, 2026).
* Added `\value` tags to `halfnormal.plot.Rd`, `plot.simplexregression.Rd`,
  and `simplexreg.methods.Rd`.
* Replaced `\dontrun{}` with `\donttest{}` in examples of `diag.im` and
  `diag.distances`; removed parallel examples that cannot run in check
  environments.
* Fixed `par()` calls in `simplexreg_plots.R` to restore graphical parameters
  via `on.exit()` immediately after modification.
* Fixed `par()` call in vignette (`relative-humidity.Rmd`) to restore
  graphical parameters after use.
* Fixed `simulate.simplexregression()` to avoid writing to `.GlobalEnv`.
* Fixed duplicated `@examples` block in `RelativeHumidity` documentation.

# SimplexRegression 0.1.1

* Reduced CRAN check time by removing heavy leave-one-out diagnostics from vignette
* Removed Wasserstein-2 distance from test suite (excessive computation time)
* Check time reduced from ~20 minutes to under 10 minutes

# SimplexRegression 0.1.0

* Initial CRAN submission.
