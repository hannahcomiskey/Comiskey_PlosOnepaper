
# Version 1.0.1
── R CMD check results ────────────────────────────────────────────────────────────────────────────────────── mcmsector 1.0.1 ────
Duration: 47.2s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔



# Version 1.0.0 ----------------------------------------------------------------
(A) Test results
Duration: 23s

❯ checking Rd line widths ... NOTE
  Rd file 'run_preprocessing_pipeline.Rd':
    \examples lines wider than 100 characters:
                                         area_classification = mcmsector::Country_and_area_classification_inclFP2020)
  
  These lines will be truncated in the PDF manual.

0 errors ✔ | 0 warnings ✔ | 1 note ✖

This NOTE is due to long example lines and does not affect functionality.
All examples run successfully.

(B) Test environments

1. macOS Sonoma 14.7.1, R 4.5.1
  R CMD check --as-cran: 
  0 errors, 0 warnings, 1 note
2. Windows (win-builder, R-release)
  Installation time in seconds: 8
  Check time in seconds: 68
  0 errors, 0 warnings, 1 note

(C) Additional comments

This is a new submission to CRAN.

(D) Notes

Possibly misspelled words in DESCRIPTION
“subnational” (22:5, 24:46)
This is not a spelling mistake. The term is intentionally used and is correct in this context.

`\examples` lines wider than 100 characters
The example line refers to stored data:
area_classification = mcmsector::Country_and_area_classification_inclFP2020

It is written on a single line for user clarity, and wrapping it would reduce readability.
This NOTE is expected and can be safely ignored.

No other issues were reported.
