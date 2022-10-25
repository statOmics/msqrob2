#' Method to fit msqrob models with robust regression and/or ridge regression and/or random effects
#' It models multiple features simultaneously, e.g. multiple peptides from the same protein.
#'
#' @description Parameter estimation of msqrob models for `QFeatures`instance.
#'              The method aggregates features within the model e.g. from peptides to proteins.
#'              It provides fold change estimates and their associated uncertainty at the aggregated
#'              level (e.g. protein level) while correcting for the peptide species that are observed
#'              in each sample. It also addresses the correlation in the data, e.g. the peptide data
#'              for the same protein in a sample are correlate because they originate from the same
#'              protein pool. The method however does not return aggregated expression values for each sample.
#'              For visualisation purposes aggregated expression values are provide by the `aggregateFeatures`
#'              function from the `QFeatures` Package
#'
#' @author Lieven Clement
#'
#' @examples
#' # Load example data
#' # The data are a Feature object with containing
#' # a SummarizedExperiment named "peptide" with MaxQuant peptide intensities
#' # The data are a subset of spike-in the human-ecoli study
#' # The variable condition in the colData of the Feature object
#' # contains information on the spike in condition a-e (from low to high)
#' data(pe)
#'
#' # Fit MSqrob model using robust ridge regression starting from peptide intensities
#' # The fold changes are calculated at the protein level while correcting for
#' # the different peptide species in each sample and the correlation between
#' # peptide intensities of peptides of the same protein in the same sample.
#' pe <- msqrobAggregate(pe, i = "peptide", fcol = "Proteins", formula = ~condition)
#' getCoef(rowData(pe[["msqrobAggregate"]])$msqrobModels[["P00956"]])
#' @param object `QFeatures` instance
#'
#' @param formula Model formula. The model is built based on the
#'     covariates in the data object.
#'
#' @param i `character` or `integer` to specify the element of the `QFeatures` that
#'        contains the log expression intensities that will be modelled.
#'
#' @param fcol The feature variable of assay ‘i’ defining how to summerise
#'        the features.
#' @param name A ‘character(1)’ naming the new assay. Default is ‘newAssay’.
#'       Note that the function will fail if there's already an assay
#'       with ‘name’.
#' @param aggregateFun A function used for quantitative feature aggregation.
#'        Details can be found in the documentation of the `aggregateFeatures`
#'        of the `QFeatures` package.
#'
#' @param modelColumnName `character` to indicate the variable name that is used
#'        to store the msqrob models in the rowData of the SummarizedExperiment
#'        instance or of the assay of the QFeatures instance. Default is "msqrobModels".
#'
#' @param robust `boolean(1)` to indicate if robust regression is
#'     performed to account for outliers. Default is `TRUE`. If
#'     `FALSE` an OLS fit is performed.
#'
#'
#' @param maxitRob `numeric(1)` indicating the maximum iterations in
#'        the IRWLS algorithm used in the M-estimation step of the robust
#'        regression.
#'
#' @param tol `numeric(1)` indicating the tolerance for declaring convergence
#'        of the M-estimation loop.
#'
#' @param doQR `boolean(1)` to indicate if QR decomposition is used when adopting
#'     ridge regression. Default is `TRUE`. If `FALSE` the predictors of the fixed
#'     effects are not transformed, and the degree of shrinkage can depend on the encoding.
#'
#' @param lmerArgs a list (of correct class, resulting from ‘lmerControl()’
#'        containing control parameters, including the nonlinear optimizer to be used
#'        and parameters to be passed through to the nonlinear optimizer, see the
#'        ‘lmerControl’ documentation of the lme4 package for more details.
#'        Default is `list(control = lmerControl(calc.derivs = FALSE))`
#'
#' @return A ‘QFeatures’ object with an additional assay.
#'
#' @rdname msqrobAggregate
#'
#' @aliases msqrobAggregate msqrobAggregate,QFeatures-method
#' @import SummarizedExperiment
#'
#' @export


setMethod(
    "msqrobAggregate", "QFeatures",
    function(object,
    formula,
    i,
    fcol,
    name = "msqrobAggregate",
    aggregateFun = MsCoreUtils::medianPolish,
    modelColumnName = "msqrobModels",
    robust = TRUE,
    ridge = FALSE,
    maxitRob = 1,
    tol = 1e-6,
    doQR = TRUE,
    lmerArgs = list(control = lmerControl(calc.derivs = FALSE))) {
        if (is.null(object[[i]])) stop("QFeatures object does not contain assay ", i)
        if (!(fcol %in% colnames(rowData(object[[i]])))) stop("The rowData of Assay ", i, " of the QFeatures object does not contain variable", fcol)
        object <- QFeatures::aggregateFeatures(
            object = object,
            i = i,
            fcol = fcol,
            name = name,
            fun = aggregateFun
        )
        
        if (length(formula) == 3) {
          formula <- formula[-2]
        }
        rowdata <- rowData(object)[[i]]
        data <- colData(object)
        
        #Get the variables from the formula and check if they are in the coldata or rowdata 
        check_vars <- all.vars(formula) %in% c(colnames(rowdata), colnames(data)) 
        if (!all(check_vars)){
          if(sum(!check_vars) >1) {
            vars_not_found <- paste0(all.vars(formula)[!check_vars], collapse=", ")
            stop(paste("Variables", vars_not_found, "are not found in coldata or rowdata"), sep = "")  
          } else{
            vars_not_found <- all.vars(formula)[!check_vars]
            stop(paste0("Variable ", vars_not_found, " is not found in coldata or rowdata"), sep = "")  
          }
        }
        
        #If there is at least one variable from the formula in the rowdata then keep rowdata
        #Otherwise return NULL
        if(any(colnames(rowdata) %in% all.vars(formula))){
          #Select only the relevant columns
          rowdata <- rowdata[colnames(rowdata) %in% all.vars(formula)]    
        } else {
          rowdata <- NULL
        }
        
        #Select only the relevant columns
        data <- data[colnames(data) %in% all.vars(formula)]
        
        rowData(object[[name]])[[modelColumnName]] <- msqrobLmer(
            y = assay(object[[i]]),
            formula = formula,
            data = colData(object),
            rowdata = rowdata,
            robust = robust,
            ridge = ridge,
            maxitRob = maxitRob,
            tol = tol,
            doQR = doQR,
            lmerArgs = lmerArgs,
            featureGroups = rowData(object[[i]])[[fcol]]
        )
        return(object)
    }
)
