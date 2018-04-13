# gpr 1.0.1

* Added a `NEWS.md` file to track changes to the package.
* Changed name of function sampleFromK to sampleFromL.
* Modified covariance function `cf` and covariance matrix function `covMat` so that noise is added when covariance matrix is generated rather than directly by the covariance function.
* Added a 0-value returning mean function `mf`
* Function `predictGP` can now take a mean function rather than assuming a 0 mean.
* Function `predictGP` now returns that original dataset.
* Function `predictGP` now returns object of class `gpr`.
* New helper function `plot.gpr` added for convenience.
* set.seed(0) added for reproducibility.
* `Vignette` and `README` modified to use new functions.
* Added `examples`
