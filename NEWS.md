# msqrob2 1.21

Changes since Bioconductor 3.23 (msqrob2 1.20).

## New features

- Empirical Bayes variance trend: argument `trend` in `msqrob()`,
  `msqrobLm()`, `msqrobLmer()`, `msqrobAggregate()` and `msqrobHurdle()`
  (intensity component). `trend = FALSE` (default) or `"none"` keeps a
  global prior variance. `"mean"` conditions the prior on mean
  log-intensity (as `limma::eBayes(trend = TRUE)`). `"count"` conditions it
  on log2 of the number of features (e.g. peptides) summarised into each
  protein: the `.n` column of `QFeatures::aggregateFeatures()` for
  protein-level data, the number of features per group for
  `msqrobAggregate()` and `msqrobLmer(featureGroups = )`, or `counts` in
  `msqrobLm()`/`msqrobLmer()`. This is similar in spirit to DEqMS. The
  binomial models (`msqrobGlm()`, `msqrobQB()`, the count component of
  `msqrobHurdle()`) have no trend.
- Designs that are rank deficient for a feature, e.g. when a condition or
  the reference level is not observed for a protein, are now fitted
  (`msqrobLm()`, `msqrobLmer()`, also with `robust = TRUE`) instead of
  returning `fitError`. Contrasts that are not estimable for that feature
  are `NA` (using the estimability package).
- `msqrobLmer()` stores the estimated variance component ratios of each fit
  in `getModel(model)$theta` (lme4's parameterisation and names).
- `msqrobLmer()` also accepts a `data.frame` as `data`.

## Changes in results

- Robust fits without random effects or ridge (`msqrob()` and the intensity
  component of `msqrobHurdle()` with such a formula, and `msqrobLm()`) use
  the variance of the Huber M-estimator (as
  `MASS::summary.rlm(method = "XtX")`) instead of the weighted least squares
  variance. The latter treats the Huber weights as known and underestimates
  the variance, which made robust tests liberal. Log2 fold changes are
  unchanged; standard errors are typically a few percent larger. Robust
  mixed and ridge models (`msqrobLmer()`, and thus `msqrobAggregate()`)
  keep the previous variance.
- Empirical Bayes variance moderation always estimates the prior with
  limma's `fitFDistUnequalDF1()` (`squeezeVar(legacy = FALSE)`), on the
  fitted features only. Previously limma switched to its moment estimator
  `fitFDist()` as soon as one feature failed to fit, so the method depended
  on fit failures and filtering. Posterior variances, and hence p-values,
  change slightly. Requires limma >= 3.62.0.
- `msqrobLm()`: `maxitRob` now defaults to 1, as in `msqrob()` and the
  other functions (was 5). Results of `msqrob()` are not affected.

## Bug fixes

- `msqrobAggregate()` (and `msqrobLmer()` with `rowdata`): when the formula
  uses exactly one sample variable and one feature variable (e.g.
  `~ condition + (1 | Sequence)`), the sample variable keeps its type (in
  the previous release a factor became integer codes and a numeric
  covariate was treated as a factor).
- `msqrobHurdle(ridge = TRUE)` stores the intensity models on the correct
  rows when the rows are not sorted by name (in the previous release it
  could fail or attach the model of another protein).
- `msqrobLmer()`: a failure after the lmer fit gives a `fitError` instead of
  an incomplete `lmer` model without coefficients.
- No error, but a warning, when no feature can be fitted.
- `msqrobLm()` no longer warns that `rlm` "failed to converge in `maxitRob`
  steps": stopping after `maxitRob` iterations is intended (a one-step
  M-estimator by default). Other warnings are kept.

## Internal

- Modular implementation of the fitting functions and unit tests for all
  functions in msqrob.R.

# msqrob2 1.20

## msqrob2 1.20.0

- Version bump 

# msqrob2 1.19

## msqrob2 1.19.1

- Add functions to calculate log-normalisation factors that can be used with 
  sweep function
- Add function msqrobCollect to collect results tables
- Add function plotVolcano to make volcanoplots for hypothesisTest tables
- Add function createPairwiseContrasts to generate contrasts for all pairwise 
  comparisons
- Illustrate these new functionalities in the vignettes
- New Vignettes for DIA (DIA-NN and Spectronaut)

## msqrob2 1.19.2

- Update cptac vignette (fread with argument check.names = TRUE)

## msqrob2 1.19.3

- Update vignettes (arguments fread check.names = TRUE and integer64 = "double") 
to avoid issues with readQFeatures
- Update README to point to msqrob2book

# msqrob2 1.17

## msqrob2 1.17.2

(nothing yet)

## msqrob2 1.17.1

- Fixed NEWS file
- Fixed vignette upon renaming of longFormat to longForm

## msqrob2 1.17.0

- New Bioconductor devel release 3.22

# msqrob2 1.16

## msqrob2 1.16.0

- New Bioconductor stable release 3.21

# msqrob2 1.15

## msqrob2 1.15.1

- fix: fixed fnames in cptac vignette and moved away from `msdata` to 
  rely on `MsDataHub` instead

## msqrob2 1.15.0

- New Bioconductor devel release 3.21

# msqrob2 1.14

## msqrob2 1.14.1

- Added more flexibility in the msqrobLM function for models with missing variables
- Fixed full rank check for models without ridge in certain cases

## msqrob2 1.14.0

- New Bioconductor stable release 3.20

# msqrob2 1.13

## msqrob2 1.13.1

- Implemented msqrobAggregate() method for SummarizedExperiment objects. 

## msqrob2 1.13.0

- New Bioconductor devel release 3.20

# msqrob2 1.12

## msqrob2 1.12.0
 
- New Bioconductor stable release 3.19

## msqrob2 1.11.2

- Fixed issue related to levels of a ridge variable
- Fixed issue when fitting only one random effect

## msqrob2 1.11.1

- Fixed issues related to reference class changes in the models
- Fixed issue related to colData assay levels
- Refactored internal code

## msqrob2 1.6.2

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
