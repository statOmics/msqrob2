
# msqrob 1.17

## msqrob 1.17.1

- Fixed NEWS file
- Fixed vignette upon renaming of longFormat to longForm

## msqrob 1.17.0

- New Bioconductor devel release 3.22

# msqrob 1.16

## msqrob 1.16.0

- New Bioconductor stable release 3.21

# msqrob 1.15

## msqrob 1.15.1

- fix: fixed fnames in cptac vignette and moved away from `msdata` to 
  rely on `MsDataHub` instead

## msqrob 1.15.0

- New Bioconductor devel release 3.21

# msqrob 1.14

## msqrob 1.14.1

- Added more flexibility in the msqrobLM function for models with missing variables
- Fixed full rank check for models without ridge in certain cases

## msqrob 1.14.0

- New Bioconductor stable release 3.20

# msqrob 1.13

## msqrob 1.13.1

- Implemented msqrobAggregate() method for SummarizedExperiment objects. 

## msqrob 1.13.0

- New Bioconductor devel release 3.20

# msqrob 1.12

## msqrob 1.12.0
 
- New Bioconductor stable release 3.19

## msqrob 1.11.2

- Fixed issue related to levels of a ridge variable
- Fixed issue when fitting only one random effect

## msqrob 1.11.1

- Fixed issues related to reference class changes in the models
- Fixed issue related to colData assay levels
- Refactored internal code

## msqrob 1.6.2

- Fixed issues related to residual degrees of freedom of overparameterized models.
- Fixed issue relating to coldata variables as factors.

## msqrob2 1.5.4

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
