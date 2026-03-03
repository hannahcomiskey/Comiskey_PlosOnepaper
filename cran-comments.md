# Version 1.0.1
> ── R CMD check results ────────────────────────────────────────────────────────────────────────────────────── mcmsector 1.0.1 ────

Duration: 47.2s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

> Test environments

* local macOS (R 4.0.0)
* Ubuntu 22.04 (R 4.0.0)
* Windows 11 (R 4.0.0)
* rhub: linux, macos, macos-arm64, windows

There was one ERROR on macOS related to JAGS.
This is due to JAGS not being available in the macOS build environment.
The package imports `rjags`, which requires a system installation of JAGS.

The package builds and checks successfully on:
- Linux 
- Windows

Since JAGS is an external system dependency, this macOS failure is not due to the package itself.

> Additional comments

This is a new submission to CRAN.
