# msqrob 1.11.1

- Fixed issues related to reference class changes in the models
- Fixed issue related to colData assay levels
- Refactored internal code

# msqrob 1.6.2

- Fixed issues related to residual degrees of freedom of overparameterized models.
- Fixed issue relating to coldata variables as factors.

# msqrob2 1.5.4

- Removed the full lmer model that was accidently left in the statmodel object.

# msqrob2 1.5.3

- Added the option to use mixed models without ridge regression
- Added the option to use rowdata variables in the mixed models

# msqrob2 1.5.1

- Fix weighted variance covariance matrix and QR decomposition in the msqrobLmer function

# msqrob2 1.1.1

- Fix filtering steps in vignette: now using `QFeatures::filterFeatures()`

# msqrob2 0.99.6

- Update authors in Description file

# msqrob2 0.99.5

- Fix standard errors on the model parameter estimates by msqrobLmer when using doQR = TRUE

# msqrob2 0.99.4

- Minor update vignette. Replace eval=FALSE in one R chunk so that the code is evaluated.

# msqrob2 0.99.3

- Resolved notes on signalers
- Resolved examples with lines wider than 100 characters
- Avoiding sapply and 1:...

# msqrob2 0.99.2

- Added a `NEWS.md` file to track changes to the package.
- Changed dependency from R (>= 4.0) to R (>= 4.1)
- Updated citation file

## Changes in version 0.99.1

 - Submission to Bioc
