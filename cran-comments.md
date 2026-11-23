## Update

This is an update from 1.0-1 to 1.1-0. It:

* fixes the NOTE on the current CRAN check results (old-style `citEntry()`
  in CITATION, now `bibentry()`),
* corrects the higher-order statistic in `hoa()` for general contrasts,
  and fixes several smaller bugs (see NEWS.md),
* moves ggplot2 from Depends to Imports, and uses suggested packages
  conditionally in examples and the vignette.

## Test environments

* local: Ubuntu 26.04, R 4.6.2
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1),
  macos-latest (R-release), windows-latest (R-release)
* win-builder: R-devel (2026-11-17 r90926 ucrt)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There is one reverse dependency, pairwiseCI. R CMD check of pairwiseCI
0.1-27 against mcprofile 1.1-0 gives no problems.
