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

## Resubmission (version 0.1.2)

This is a resubmission addressing the comments from Konstanze Lauseker.
All requested changes have been made:

* Removed redundant "Provides functions for" from DESCRIPTION.
* Added references (Barndorff-Nielsen & Jorgensen, 1991; Justino &
  Cribari-Neto, 2026) to DESCRIPTION in the required format.
* Added \value tags to halfnormal.plot.Rd, plot.simplexregression.Rd,
  and simplexreg.methods.Rd.
* Replaced \dontrun{} with \donttest{} in examples of diag.im and
  diag.distances; parallel examples were removed as they cannot run
  in check environments.
* Fixed par() calls in R/simplexreg_plots.R to use on.exit() immediately
  after modification.
* Fixed par() call in vignette to restore graphical parameters after use.
* Fixed simulate.simplexregression() to avoid writing to .GlobalEnv.
