#' Function to fit msqrob models using lm and rlm
#'
#' @description Low-level function for parameter estimation with msqrob
#'              using the ordinary least squares or robust regression
#'              base on the MASS::rlm function.
#'
#' @param y A `matrix` with the quantified feature intensities. The
#'        features are along the rows and samples along the columns.
#'
#' @param formula Model formula. The model is built based on the
#'        covariates in the data object.
#'
#' @param data A `DataFrame` with information on the design. It has
#'        the same number of rows as the number of columns (samples) of
#'        `y`.
#'
#' @param robust `boolean(1)` to indicate if robust regression is
#'        performed to account for outliers. Default is `TRUE`. If
#'        `FALSE` an OLS fit is performed.
#'
#' @param maxitRob `numeric(1)` indicating the maximum iterations in
#'        the IRWLS algorithm used in the M-estimation step of the robust
#'        regression.
#'
#' @examples
#'
#' # Load example data
#' # The data are a Feature object with containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#'
#' # Aggregate peptide intensities in protein expression values
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#' pe
#'
#' # Fit MSqrob model using robust regression with the MASS rlm function
#' models <- msqrobLm(assay(pe[["protein"]]), ~condition, colData(pe))
#' #' getCoef(models[[1]])
#' @return A list of objects of the `StatModel` class.
#'
#' @rdname msqrobLm
#'
#' @author Lieven Clement, Oliver M. Crook
#'
#' @importFrom MASS rlm
#' @importFrom stats model.matrix lm.fit
#' @importFrom limma squeezeVar
#' @importFrom methods is
#'
#' @export
msqrobLm <- function(y,
    formula,
    data,
    robust = TRUE,
    maxitRob = 5) {
    data <- .matchQuantColsOrder(data, y)
    myDesign <- model.matrix(formula, data)
    # apply the model to each protein
    models <- apply(y, 1,
        function(y, design) {
            ## computability check
            obs <- is.finite(y)
            type <- "fitError"
            model <- list(
                coefficients = NA, vcovUnscaled = NA,
                sigma = NA, df.residual = NA, w = NA
            )

            if (sum(obs) > 0) {
                ## subset to finite observations, attention with R column switching
                X <- design[obs, , drop = FALSE]
                X <- X[,colMeans(X == 0) != 1 , drop = FALSE]
                y <- y[obs]
                colnames_orig <- colnames(design)

                if (robust) {
                    ## use robust regression from MASS package, "M" estimation is used
                    mod <- try(MASS::rlm(X, y,
                        method = "M",
                        maxit = maxitRob
                    ),
                    silent = TRUE
                    )
                    if (!is(mod, "try-error")) {
                        type <- "rlm"
                    }
                } else {
                    ## if robust regression is not performed use standard linear fit
                    mod <- try(lm.fit(X, y))
                    if ((!is(mod, "try-error")) & mod$rank == ncol(X)) {
                        type <- "lm"
                    }
                }

                if (type == "rlm") {
                    w <- mod$w
                    sigma <- sqrt(sum(mod$w * mod$resid^2) / (sum(mod$w) - mod$rank))
                    df.residual <- sum(mod$w) - mod$rank
                    if (df.residual < 2L) type <- "fitError"
                }

                if (type == "lm") {
                    w <- NULL
                    sigma <- sqrt(sum(mod$residuals^2 / mod$df.residual))
                    df.residual <- mod$df.residual
                    if (df.residual < 2L) type <- "fitError"
                }

                if (type != "fitError") {
                    coef <- rep(NA, length(colnames_orig))
                    names(coef) <- colnames_orig
                    coef[names(mod$coef)] <- mod$coef
                    vcovUnscaled <- matrix(NA, nrow =length(colnames_orig), ncol = length(colnames_orig))
                    rownames(vcovUnscaled) <- colnames(vcovUnscaled) <-  colnames_orig
                    vcovUnscaled[names(mod$coef), names(mod$coef)] <- msqrob2:::.vcovUnscaled(mod)

                    model <- list(
                        coefficients = coef,
                        vcovUnscaled = .vcovUnscaled(mod),
                        sigma = sigma,
                        df.residual = df.residual,
                        w = w
                    )
                }
            }
            ## return object of class Statmodel (from apply)
            .StatModel(
                type = type,
                params = model,
                varPosterior = as.numeric(NA),
                dfPosterior = as.numeric(NA)
            )
        },
        design = myDesign
    ) ## end of apply here

    ## Squeeze a set of sample variances together by computing
    ## empirical Bayes posterior means
    hlp <- limma::squeezeVar(
        var = vapply(models, getVar, numeric(1)),
        df = vapply(models, getDF, numeric(1))
    )

    ## Put variance and degrees of freedom in appropriate slots
    for (i in seq_len(length(models))) {
        mydf <- hlp$df.prior + getDF(models[[i]])
        models[[i]]@varPosterior <- as.numeric(hlp$var.post[[i]])
        models[[i]]@dfPosterior <- as.numeric(mydf)
    }

    ## Return object of class StatModel
    return(models)
}


#' Function to fit msqrob models with ridge regression and/or random effects using lme4
#'
#' @description Low-level function for parameter estimation with msqrob
#'              using the robust ridge regression. The models can be fitted for each
#'              feature (e.g. summarised protein expression values) or multiple features
#'              belonging to the same accession can be modelled simultaneously
#'              e.g. peptide-based models where all peptide intensities for the same
#'              protein are modelled simultaneously. The fold changes and uncertainty
#'              estimates are then calculated at the protein level while correcting
#'              for peptide species and within sample correlation.
#'
#' @param y A `matrix` with the quantified feature intensities. The
#'        features are along the rows and samples along the columns.
#'
#' @param formula Model formula. The model is built based on the
#'        covariates in the data object.
#'
#' @param data A `DataFrame` with information on the design. It has
#'        the same number of rows as the number of columns (samples) of
#'        `y`.
#'
#' @param rowdata A `DataFrame` with the rowData information of the SummarizedExperiment.
#'        It has the same number of rows as the number of rows (features) of
#'        `y`.
#'
#' @param robust `boolean(1)` to indicate if robust regression is
#'        performed to account for outliers. Default is `TRUE`. If
#'        `FALSE` an OLS fit is performed.
#'
#' @param ridge `boolean(1)` to indicate if ridge regression is
#'        performed. Default is `FALSE`. If `TRUE` the fixed effects are
#'        estimated via penalized regression and shrunken to zero.
#'
#' @param maxitRob `numeric(1)` indicating the maximum iterations in
#'        the IRWLS algorithm used in the M-estimation step of the robust
#'        regression.
#'
#' @param tol `numeric(1)` indicating the tolerance for declaring convergence
#'        of the M-estimation loop.
#'
#'
#' @param doQR `boolean(1)` to indicate if QR decomposition is used when adopting
#'        ridge regression. Default is `TRUE`. If `FALSE` the predictors of the fixed
#'        effects are not transformed, and the degree of shrinkage can depend on the encoding.
#'
#' @param featureGroups vector of type `character` or vector of type `factor` indicating how to aggregate
#'        the features. Is only used when multiple features are used to build the model, e.g. when starting
#'        from peptide data and modelling the fold change at the protein level. The default is `NULL`
#'
#' @param lmerArgs a list (of correct class, resulting from ‘lmerControl()’
#'        containing control parameters, including the nonlinear optimizer to be used
#'        and parameters to be passed through to the nonlinear optimizer, see the
#'        ‘lmerControl’ documentation of the lme4 package for more details.
#'        Default is `list(control = lmerControl(calc.derivs = FALSE))`
#'
#' @examples
#'
#' # Load example data
#' # The data are a Feature object with containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#'
#' # Aggregate peptide intensities in protein expression values
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#'
#' # Fit MSqrob model using robust ridge regression upon summarization of
#' # peptide intensities into protein expression values
#' modelsRidge <- msqrobLmer(assay(pe[["protein"]]), ~condition, data = colData(pe),
#'                           ridge = TRUE)
#' getCoef(modelsRidge[[1]])
#'
#' # Fit MSqrob model using robust ridge regression starting from peptide intensities
#' # The fold changes are calculated at the protein level while correcting for
#' # the different peptide species in each sample and the correlation between
#' # peptide intensities of peptides of the same protein in the same sample.
#' # Add the samples variable to colData
#' colData(pe)$samples <- rownames(colData(pe))
#' modelsPepBased <- msqrobLmer(assay(pe[["peptide"]]),
#'     formula = ~condition + (1|samples) + (1|Sequence), data = colData(pe),
#'     rowdata = rowData(pe[["peptide"]]), featureGroups = rowData(pe[["peptide"]])$Proteins,
#'     ridge = TRUE)
#' getCoef(modelsPepBased[[1]])
#' @return A list of objects of the `StatModel` class.
#'
#' @rdname msqrobLmer
#'
#' @author Lieven Clement, Oliver M. Crook
#'
#' @importFrom MASS psi.huber
#' @importFrom stats resid update.formula resid mad
#' @importFrom methods as is
#' @import lme4
#' @import Matrix
#' @importFrom BiocParallel bplapply bpmapply
#' @importFrom MultiAssayExperiment DataFrame
#'
#' @export

msqrobLmer <- function(y,
                       formula,
                       data,
                       rowdata = NULL,
                       tol  =1e-6,
                       robust = TRUE,
                       ridge = FALSE,
                       maxitRob = 1,
                       doQR = TRUE,
                       featureGroups=NULL,
                       lmerArgs = list(control = lmerControl(calc.derivs = FALSE)),
                       keep.model = FALSE) {

  #Get the featureGroups variable
  if (is.null(featureGroups)){
    featureGroups <- rownames(y)
  }

  if (!is.null(rowdata)){
    #select only the relevant columns
    rowdata <- rowdata[colnames(rowdata) %in% all.vars(formula)]
    rowdata <- split.data.frame(rowdata, featureGroups)
  }

  #Select only the relevant columns
  data <- data[,colnames(data) %in% all.vars(formula), drop = FALSE]

  y <- split.data.frame(y, featureGroups)

  fit_args <- list(
    formula = formula,
    coldata = data,
    robust = robust,
    maxitRob = maxitRob,
    tol = tol,
    keep.model = keep.model
  )
  if (ridge) {
    fit_func <- .ridge_msqrobLmer
    fit_args$doQR <- doQR
  } else {
    fit_func <- .noridge_msqrobLmer
  }
  models <- if(is.null(rowdata)){
    do.call(bplapply, append(list(y, FUN=fit_func), fit_args))
  } else {
    bpmapply(FUN = fit_func, y, rowdata, MoreArgs = fit_args)
  }

  modelVars <- vapply(models, getVar, numeric(1))
  modelDfs <- vapply(models, getDF, numeric(1))
  modelMask <- !is.na(modelVars) & !is.na(modelDfs)
  varPosterior <- rep_along(models, NA_real_)
  hlp <- limma::squeezeVar(
    var = modelVars[modelMask],
    df = modelDfs[modelMask],
  )
  varPosterior[modelMask] <- hlp$var.post

  for (i in seq_len(length(models))) {
    models[[i]]@varPosterior <- as.numeric(varPosterior[[i]])
    models[[i]]@dfPosterior <- as.numeric(hlp$df.prior + getDF(models[[i]]))
  }
  return(models)
}

## Fit the mixed models with ridge regression
.ridge_msqrobLmer <- function(y, rowdata=NULL, formula, coldata, doQR,
                              robust, maxitRob=1, tol = 1e-06,
                              keep.model = FALSE){

  #Create the matrix containing the variable information
  data <- .create_data(y,rowdata,coldata)

  #all necessary variables are now in data,  now we can create the fixed object if we use ridge regression
  fixed <- model.matrix(nobars(formula), data = data)
  data$fixed <- fixed
  data$y <- as.matrix(y)
  data <- data[!is.na(data$y), , drop = FALSE]

  #Checking reference class changes
  #nonestimable_paramaters <- limma::nonEstimable(data$fixed)
  data$fixed <- data$fixed[,colMeans(data$fixed == 0) != 1 , drop = FALSE]

  if (sum(!grepl("(Intercept)", colnames(fixed))) < 2 & nobars(formula)[[2]] != 1) {
    stop("The mean model must have more than two parameters for ridge regression.
              if you really want to adopt ridge regression when your factor has only two levels
              rerun the function with a formula where you drop the intercept. e.g. ~-1+condition
            ")
  }

  if(is.null(findbars(formula))) {
    formula <- formula(y ~ (1|ridge))
  } else {
    if (nobars(formula)[[2]] != ~1){
      # update formula to remove any fixed effect variables and replace with ridge
      formula <- formula(
        paste0("y ~ (1|ridge) + ", paste0("(",paste(findbars(formula), collapse=")+("),")")))
    } else {
      formula <- update.formula(formula, y~.)
    }
  }

  qrFixed <- qr(data$fixed)

  if (doQR) {
    Q <- qr.Q(qrFixed)
  } else {
    Q <- data$fixed
  }

  model <- NULL
  ##Fooling lmer to adopt ridge regression using Fabian Scheipl's trick
  if (qrFixed$rank == ncol(data$fixed)){
    try({
      colnames(Q) <- colnames(data$fixed)
      if (colnames(data$fixed)[1] == "(Intercept)") {
        Q <- Q[, -1]
      }

      data$ridge <- factor(rep(colnames(Q), length = nrow(data)), levels = colnames(Q))

      #Parse the data and formula
      parsedFormulaC <- lFormula(formula,data = as.list(data))
      parsedFormulaC$reTrms$cnms$ridge <- ""
      ridgeId <- grep(names(parsedFormulaC$reTrms$Ztlist), pattern = "ridge")
      parsedFormulaC$reTrms$Ztlist[[ridgeId]] <- as(Matrix(t(Q)), class(parsedFormulaC$reTrms$Ztlist[[ridgeId]]))
      parsedFormulaC$reTrms$Zt <- do.call(rbind,parsedFormulaC$reTrms$Ztlist)

      #Create deviance function to be optimized
      devianceFunctionC <- do.call(mkLmerDevfun, parsedFormulaC)
      #optimize deviance function
      optimizerOutputC <- optimizeLmer(devianceFunctionC)
      #Package up the results
      model <- mkMerMod(rho = environment(devianceFunctionC),
                        opt = optimizerOutputC,
                        reTrms = parsedFormulaC$reTrms,
                        fr = parsedFormulaC$fr)
    }, silent=TRUE)
  }

  if (is.null(model)) {
    type <- "fitError"
    model <- list(coefficients = NA, vcovUnscaled = NA, sigma = NA, df.residual = NA)
  } else {
    df.residual <- 0
    try({
      type <- "lmer"
      #extract deviance composition
      #pwrss=penalied weighted residual sum of squares
      model@frame$`(weights)` <- rep(1, dim(model@frame)[1])
      sseOld <- model@devcomp$cmp['pwrss']
      if (robust == TRUE){
        model <- .robust_fitting(model, maxitRob, sseOld, tol)
      }

      sigma <- sigma(model)
      betas <- .getBetaB(model)
      vcovUnscaled <- as.matrix(.getVcovBetaBUnscaled(model))
      if (nobars(formula)[[2]] != 1) {
        if (doQR) {
          if (colnames(data$fixed)[1] == "(Intercept)") {

            ids <- c(1, grep("ridge", names(betas)))
          } else {

            ids <- grep("ridge", names(betas))
          }

          Rinv <- diag(length(betas))
          coefNames <- names(betas)
          Rinv[ids, ids] <- solve(qr.R(qrFixed))
          Rinv[1, 1] <- 1
          betas <- c(Rinv %*% betas)
          names(betas) <- coefNames

          vcovUnscaled <- Rinv %*% vcovUnscaled %*% t(Rinv)
          rownames(vcovUnscaled) <- colnames(vcovUnscaled) <- names(betas)
        }
      }

      df.residual <- .getDfLmer(model)
      if(is.na(df.residual)){
        df.residual <- 0
      }
    }, silent = TRUE)

    model <- .create_model(betas, vcovUnscaled, sigma, df.residual, w, model, keep.model = keep.model)
  }

  return(StatModel(type = type,
                   params = model,
                   varPosterior = NA_real_,
                   dfPosterior = NA_real_))
}

## Fit the mixed models without ridge regression
.noridge_msqrobLmer <- function(y, rowdata=NULL, formula, coldata,
                                robust, maxitRob=0, tol = 1e-06,
                                keep.model = FALSE){
  #Create the matrix containing the variable information
  data <- .create_data(y,rowdata,coldata)

  data_model_matrix <- model.matrix(nobars(formula), data = data)
  formula <- update.formula(formula, y~.)

  data$y <- as.matrix(y)
  data_model_matrix <- data_model_matrix[!is.na(data$y), , drop = FALSE]
  data <- data[!is.na(data$y), , drop = FALSE]
  #Checking for reference class changes
  #nonestimable_parameters <- limma::nonEstimable(data_model_matrix)
  data_model_matrix <- data_model_matrix[,colMeans(data_model_matrix == 0) != 1 , drop = FALSE]
  model <- NULL

  if(qr(data_model_matrix)$rank == ncol(data_model_matrix)){
    try({
      model <- lmer(formula,  as.data.frame(data))
    }, silent=TRUE)
  }

  if (is.null(model)) {
    type <- "fitError"
    model <- list(coefficients = NA, vcovUnscaled = NA, sigma = NA, df.residual = NA)
  } else {
    df.residual <- 0
    try({
      type <- "lmer"
      #extract deviance composition
      #pwrss=penalied weighted residual sum of squares
      model@frame$`(weights)` <- rep(1, dim(model@frame)[1])
      sseOld <- model@devcomp$cmp['pwrss']

      if (robust == TRUE){
        model <- .robust_fitting(model, maxitRob, sseOld, tol)
      }

      sigma <- sigma(model)
      betas <- .getBetaB(model)
      vcovUnscaled <- .getVcovBetaBUnscaled(model)
      df.residual <- .getDfLmer(model)
      if(is.na(df.residual)){
        df.residual <- 0
      }
    }, silent = TRUE)

    model <- .create_model(betas, vcovUnscaled, sigma, df.residual, w,
                           model, keep.model = keep.model)
  }

  return(StatModel(type = type,
                   params = model,
                   varPosterior = NA_real_,
                   dfPosterior = NA_real_))
}


## Calculate unscaled covariance matrix for lm or rlm fit
.vcovUnscaled <- function(model) {
    p1 <- 1L:model$rank
    p <- length(model$coefficients)

    out <- matrix(NA, p, p)
    out[!is.na(model$coefficients), !is.na(model$coefficients)] <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
    colnames(out) <- rownames(out) <- names(model$coefficients)

    return(out)
}

#' @import lme4
#' @import Matrix
#' @importFrom methods cbind2


.getVcovBetaBUnscaled <- function(model) {
    X <- lme4::getME(model, "X")
    Z <- lme4::getME(model, "Z")
    XZ <- cbind2(X,Z)

    if (is.null(model@frame$`(weights)`)){
      model@frame$`(weights)` <- 1
    }

    vcovInv <- Matrix::crossprod(model@frame$`(weights)`^.5 * XZ)
    Ginv <- Matrix::solve(
        Matrix::tcrossprod(getME(model, "Lambda")) +
            Matrix::Diagonal(ncol(Z), 1e-18)
    )

    i <- -seq_len(ncol(X))
    vcovInv[i, i] <- vcovInv[i, i] + Ginv
    vcovInv <- Matrix::solve(vcovInv)
    ranefLevels <- imap(model@flist, ~ {
        paste0(.y, levels(.x))
    })
    zNames <- unlist(lapply(seq_len(length(model@cnms)),
        function(x, cnms, levels) {
            c(outer(cnms[[x]], levels[[names(cnms)[x]]], paste0))
        },
        cnms = model@cnms,
        levels = ranefLevels
    ))

    rownames(vcovInv) <- colnames(vcovInv) <- c(colnames(X), zNames)
    return(vcovInv)
}

#' @import purrr
#' @import lme4
.getBetaB <- function(model) {
    betaB <- c(as.vector(lme4::getME(model, "beta")), as.vector(lme4::getME(model, "b")))
    ranefLevels <- purrr::imap(model@flist, ~ {
        paste0(.y, levels(.x))
    })
    zNames <- unlist(lapply(seq_len(length(model@cnms)), function(x, cnms, levels) {
        c(outer(cnms[[x]], levels[[names(cnms)[x]]], paste0))
    },
    cnms = model@cnms, levels = ranefLevels
    ))
    names(betaB) <- c(colnames(model@pp$X), zNames)
    betaB
}

#' Calculate the weighted REML residual degrees of freedom
#' @importFrom stats resid
.getDfLmer <- function(object) {
    w <- object@frame$"(weights)"
    if (is.null(w)) w <- 1
    sigma <- sigma(object)
    sum((resid(object) * sqrt(w))^2) / sigma^2
}

.robust_fitting <- function(model, maxitRob, sseOld, tol){
  while (maxitRob > 0) {
    maxitRob <- maxitRob - 1
    res <- resid(model)
    model@frame$`(weights)` <- MASS::psi.huber(res / (mad(res, 0)))
    model <- refit(model)
    sse <- model@devcomp$cmp["pwrss"]
    if (abs(sseOld - sse) / sseOld <= tol) break
    sseOld <- sse
  }
  return(model)
}

.create_model <- function(betas, vcovUnscaled, sigma, df.residual, w, model, keep.model = FALSE){
  if (df.residual<2L){
    res <- list(coefficients = NA,
                vcovUnscaled = NA,
                sigma = NA,
                df.residual = NA,
                w = NA)
  } else {
    res <- list(coefficients = betas,
                vcovUnscaled = vcovUnscaled,
                sigma = sigma,
                df.residual = df.residual,
                w = model@frame$`(weights)`)
  }
  if (keep.model) {
    res$model <- model
  }
  return(res)
}


.create_data <- function(y,rowdata,coldata){
  if (is.null(rowdata)){
    data <- coldata[rep(1:nrow(coldata), each = nrow(y)), , drop = FALSE]
  } else {
    data <- cbind(
      coldata[rep(1:nrow(coldata), each = nrow(y)), ],
      rowdata[rep(1:nrow(rowdata), ncol(y)),]
    )
    data <- DataFrame(data)
    colnames(data) <- c(colnames(coldata),colnames(rowdata))
  }
  return(data)
}

#' Function to fit msqrob models to peptide counts using glm
#'
#' @description Low-level function for parameter estimation with msqrob
#'              by modeling peptide counts using quasibinomial glm
#'
#' @param y A `matrix` with the peptide counts. The
#'        features are along the rows and samples along the columns.
#'
#' @param npep A vector with number of peptides per protein. It has as length
#'        the number of rows of y. The counts are equal or larger than the largest
#'        peptide count in y.
#'
#' @param formula Model formula. The model is built based on the
#'        covariates in the data object.
#'
#' @param data A `DataFrame` with information on the design. It has
#'        the same number of rows as the number of columns (samples) of
#'        `y`.
#'
#' @param priorCount A 'numeric(1)', which is a prior count to be added to the observations to shrink
#'          the estimated log-fold-changes towards zero.
#'
#' @param binomialBound logical, if ‘TRUE’ then the quasibinomial variance estimator will
#'        be never smaller than 1 (no underdispersion).
#'
#' @examples
#'
#' # Load example data
#' # The data are a Feature object with containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#'
#' # Aggregate peptide intensities in protein expression values
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#' pe
#'
#' # Fit MSqrob model using robust regression with the MASS rlm function
#' models <- msqrobGlm(
#'     aggcounts(pe[["protein"]]),
#'     rowData(pe[["protein"]])[[".n"]],
#'     ~condition,
#'     colData(pe)
#' )
#' getCoef(models[[1]])
#' @return A list of objects of the `StatModel` class.
#'
#' @rdname msqrobGlm
#'
#' @author Lieven Clement
#'
#' @importFrom limma squeezeVar
#' @importFrom stats model.matrix glm.fit binomial
#' @importFrom methods is
#'
#' @export

msqrobGlm <- function(y,
    npep,
    formula,
    data,
    priorCount = .1,
    binomialBound = TRUE) {
    myDesign <- model.matrix(formula, data)
    models <- lapply(seq_len(nrow(y)),
        function(i, y, npep, myDesign) {
            type <- "fitError"
            model <- list(
                coefficients = NA, vcovUnscaled = NA,
                sigma = NA, df.residual = NA, w = NULL
            )
            if (npep[i] >= max(y[i, ])) {
                mod <- try(glm.fit(
                    y = cbind(y[i, ], npep[i] - y[i, ]) + priorCount,
                    x = myDesign,
                    family = binomial()
                ))
                if ((!is(mod, "try-error")) & mod$rank == ncol(myDesign)) {
                    type <- "quasibinomial"
                }
            }
            if (!is(mod, "try-error")) {
                if (mod$deviance < 0) mod$deviance <- sum(pmax(mod$family$dev.resids(mod$y, mod$fitted.values, mod$prior.weights), 0))
                if (mod$df.residual < 2L) type <- "fitError"
            }
            if (type != "fitError") {
                model <- list(
                    coefficients = mod$coef,
                    vcovUnscaled = .vcovUnscaled(mod),
                    sigma = sqrt(mod$deviance / mod$df.residual),
                    df.residual = mod$df.residual,
                    w = mod$w
                )
            }

            ## return object of class Statmodel (from apply)
            .StatModel(
                type = type,
                params = model,
                varPosterior = as.numeric(NA),
                dfPosterior = as.numeric(NA)
            )
        },
        y = y, npep = npep, myDesign = myDesign
    )
    hlp <- limma::squeezeVar(
        var = vapply(models, getVar, numeric(1)),
        df = vapply(models, getDF, numeric(1))
    )

    for (i in seq_len(length(models))) {
        models[[i]]@varPosterior <- as.numeric(hlp$var.post[i])
        models[[i]]@dfPosterior <- as.numeric(hlp$df.prior + getDF(models[[i]]))

        if (!is.na(models[[i]]@varPosterior) & binomialBound) {
            if (models[[i]]@varPosterior < 1) {
                models[[i]]@varPosterior <- 1
                models[[i]]@dfPosterior <- Inf
            }
        }
    }

    return(models)
}
