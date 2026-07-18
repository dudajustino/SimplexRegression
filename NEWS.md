# SimplexRegression 0.1.5

## New features
* Added the `AISRowing` dataset.
* Added `r2()`, a function to compute [pseudo-]R² for `simplexregression`
  model fits.

## Bug fixes, validation, and breaking changes
* Added validation of `lambda > 0` in `parametric_mean_links`.
* Added validation of `kappa > 0` in `penalized.ss()` and `penalized.ic()`.
* Added error handling in functions that use `solve()`.
* Removed internal `par()` calls from `plot.simplexregression()`, 
  `halfnormal.plot()`, `local.influence()`, `diag.im()` and `diag.distances()`;
  these functions no longer override the user's  graphical parameters. 
  As a result, the `reset.par` argument in `plot.simplexregression()` was removed, 
  as it is no longer needed.
  
## Documentation
* Standardized and shortened variable names in the `Biomass`,
  `RelativeHumidity`, and `AbortionOpposition` datasets for consistency
  across package datasets.

## Dependencies
* Added `moments` and `parallel` to Imports in DESCRIPTION.

# SimplexRegression 0.1.4
 
* Fixed `opt$counts` handling in `simplexreg.fit()`.
* Added a `digits` argument to `penalized.ss()` and `penalized.ic()`.
* `plot.simplexregression()` now accepts a cutoff/threshold argument for
  flagging observations in residual plots.
* Added the `AbortionOpposition` dataset.
* Standardized variable names across package datasets.

# SimplexRegression 0.1.3

* Initial CRAN release.
