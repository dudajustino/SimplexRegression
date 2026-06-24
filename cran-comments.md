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

## Resubmission (version 0.1.3)

This is a resubmission addressing the remaining comments from 
Konstanze Lauseker regarding graphical parameter restoration.

All requested changes have been made:

* Fixed `par()` call in `plot.simplexregression.Rd` example to properly
  save and restore graphical parameters using `oldpar <- par(...)` and
  `par(oldpar)`.
* Fixed `options()` call in vignette (`relative-humidity.Rmd`) to properly
  save and restore user settings using `old_opts <- options(...)` and
  `options(old_opts)`.
* Recompiled vignettes to ensure both `par()` and `options()` fixes are
  reflected in `doc/relative-humidity.R`.
* Re-roxygenized documentation to ensure all `.Rd` changes are reflected.
* All graphical parameters modified in examples, vignettes, and demos
  are now properly restored to user's original settings.
