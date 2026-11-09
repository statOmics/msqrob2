# msqrob 1.21 
## msqrob 1.21.2

- `msqrobLmer()` stores the estimated variance component ratios of each
  lmer fit in `getModel(model)$theta` (lme4's parameterisation and names).
- Failed fits (`fitError`) get `NA` instead of a log2 fold change of 0 in
  `hypothesisTest()` and `topFeatures()`.
- `msqrobLmer()` no longer crashes on a feature without any observations; such
  features become `fitError`.
- `msqrobLmer()` also accepts a `data.frame` as `data`.
- Ridge models: the mean-model parameters keep the `ridge` prefix again (e.g.
  `ridgeconditionB`), so contrasts from
  `createPairwiseContrasts(ridge = TRUE)` and existing ridge contrasts work
  again (1.21.1 gave a log2 fold change of 0).
- `getContrast()` and `varContrast()` give `NA` instead of 0 for contrasts on
  names that are no parameters of the model (e.g. a typo).
- `msqrobLmer()` fits mixed models with a rank-deficient design (e.g. a
  missing reference level) instead of returning `fitError`; non-estimable
  contrasts are `NA`, as for `msqrobLm()`.
- `msqrobAggregate()` (and `msqrobLmer()` with `rowdata`): when the formula
  uses exactly one sample variable and one feature variable (e.g.
  `~ condition + (1 | Sequence)`), the sample variable lost its type: a factor
  became integer codes and a numeric covariate was treated as a factor.
- `msqrobHurdle(ridge = TRUE)`: intensity models are stored on the correct
  rows (when the rows were not sorted by name, proteins got the intensity
  model of another protein).
- `trend = "count"` with random effects or ridge: the counts are linked to the
  correct features; missing counts give an error before fitting.
- `trend = "mean"` no longer crashes when a feature has no observations.
- `trend = TRUE` gives an informative error (use `trend = "mean"`).
- `msqrobLmer()`: a failure after the lmer fit gives a `fitError`, no longer an
  incomplete `lmer` model.
- No error, but a warning, when no feature can be fitted; failed fits
  (`fitError`) keep `NA` posterior variances and degrees of freedom.
- Empirical Bayes variance moderation now always estimates the prior with
  limma's `fitFDistUnequalDF1()` (`squeezeVar(legacy = FALSE)`), using the
  fitted features only. Previously, limma used its moment estimator
  `fitFDist()` as soon as one feature failed to fit, and `fitFDistUnequalDF1()`
  otherwise, so the method depended on fit failures and filtering; the moment
  estimator can also give an infinite prior df when residual df are small.
  Posterior variances, and hence p-values, change slightly. Requires
  limma >= 3.62.0.

## msqrob 1.21.1 

- Modular implementation of msqrob.R 
- Add estimability functionality to catch and solve issues with changing 
  reference class 
- Adding unit tests for all functions in msqrob.R
- Trended squeezevar using mean (limma trend) or using feature counts used for 
  aggregation (DEqMS). 

## msqrob 1.21.0 

- Version bump

# msqrob 1.20 

## msqrob 1.20.0 

- Version bump 

# msqrob 1.19

## msqrob 1.19.1

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

# msqrob 1.17

## msqrob 1.17.2

(nothing yet)

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
