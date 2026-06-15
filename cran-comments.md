## Test environments
* local Windows 11, R 4.5.0
* win-builder (R-devel)
* GitHub Actions (ubuntu-latest, macos-13, windows-latest) with R-devel

## R CMD check results
0 errors | 0 warnings | 1 note

* This is a new submission.

* Possibly misspelled words in DESCRIPTION: "simplex", "Scout", 
  "submodel", "Cribari-Neto", "Cook's" — these are correct 
  statistical terms and proper names.

## Reverse dependencies
There are no reverse dependencies.

## Resubmission information (version 0.1.1)

This is a resubmission of SimplexRegression 0.1.0.

### Changes made to address CRAN check time NOTE:

#### Test suite modifications:
* Removed Wasserstein-2 (W2) distance from tests — this distance required
  numerical integration (quadgk) and leave-one-out refits for each
  observation, causing excessive computation time

#### Vignette modifications:
* Removed heavy diagnostic calls (`diag.im` and `diag.distances`) from the
  vignette because they perform leave-one-out refits (n = 312 observations)
  with numerical integration, which previously caused check times > 20 minutes
* All model fitting, summaries, hypothesis tests, and lightweight diagnostics
  are preserved; the vignette remains fully illustrative
* The removed functions are still available in the package for users who
  wish to apply them locally

### Impact on check time:
* Previous check time: ~20 minutes (exceeded CRAN's 10-minute limit)
* Current check time: under 10 minutes (tested locally)

## Additional comments
This package implements simplex regression models with flexible link functions.
All tests pass and documentation is complete. The package has been tested on
Windows, macOS, and Linux platforms.
