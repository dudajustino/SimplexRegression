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

## Additional notes
This submission adds `moments` and `parallel` to Imports. It also adds
input validation (lambda > 0, kappa > 0) and error handling around
solve() calls to prevent uninformative failures on invalid input or
singular matrices, in addition to new features and documentation
improvements. No previous CRAN comments are being addressed in this
submission.
