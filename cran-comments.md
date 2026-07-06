## Test environments
* local Windows 11, R 4.5.0
* win-builder (R-devel)
* GitHub Actions (ubuntu-latest, macos-13, windows-latest) with R-devel

## R CMD check results
0 errors | 0 warnings | 1 note

* Possibly misspelled words in DESCRIPTION: "simplex", "Scout", 
  "submodel", "Cribari-Neto", "Cook's" — these are correct 
  statistical terms and proper names.

## Reverse dependencies
There are no reverse dependencies.

## Summary of changes (version 0.1.4)
 
This is a routine update with bug fixes, new features, and documentation
improvements. No previous CRAN comments are being addressed in this
submission.
 
* Fixed `opt$counts` handling in `simplexreg.fit()`.
* Made decimal-place formatting consistent across `summary()`, `print()`,
  and `print.summary()` methods.
* Added a `digits` argument to `penalized.ss()` and `penalized.ic()`.
* `plot.simplexregression()` now accepts a cutoff/threshold argument for
  flagging observations in residual plots.
* Improved English wording and mathematical notation in function
  documentation.
* Added a new dataset, `AbortionOpposition`.
* Standardized variable names across package datasets (converted to
  lowercase).
