#' Function to fit msqrob models
#'
#' @description Parameter estimation of msqrob models.
#'
#' @param pe Proteomics Experiment Object, is an object off class Feature set that contains the assay with the quantified MS intensities.
#' @param assay integer or string referring to the assay of the Feature object containing the intensities for modelling, e.g. peptide or protein.
#' @param formula Model formula. The model is built based on the covariates in the colData of the pe object.
#' @param robust bolean to indicate if robust regression is performed to account for outliers. If 'FALSE' an OLS fit is performed.
#' @param maxitRob Maximum iterations in the IRWLS algorithm used in the M-estimation step of the robust regression.
#' @examples
#' @return A list of objects of the StatModel class
#' @rdname msqrob
#' @author Lieven Clement, Oliver M. Crook
#' @export
msqrob <- function(pe, 
                   assay, 
                   formula, 
                   robust = TRUE, 
                   maxitRob = 5)
{
# extract quantitative assay data  
yAll <- assay(pe[[assay]])
ngenes <- nrow(yAll) # the number of features
myDesign <- model.matrix(formula, colData(pe))

models <- apply(yAll, 1, function(y, design)
  {
    # computatability ceck
    obs <- is.finite(y)
    type <- "fitError"
    model <- list()
    if (sum(obs) > 0) {
    
     # subset to finite observations, attention with R column switching
     X <- design[obs, , drop = FALSE]
     y <- y[obs]
     
     if(robust) {
    
     # use robust regresssion from MASS package, "M" Setimation is used
     model <- try(MASS::rlm(X,
                            y, 
                            method = "M",
                            maxit = maxitRob),
                            silent = TRUE)
      if (class(model)[1]=="try-error"){ 
        
        # catch error
        model <- list()
      } else {
        
        type <- "rlm"
        class(model) <- "list"
      }
     } else {
     # if robust regression is not performed use standard linear fit   
     model <- lm.fit(X, y)
     type <- "lm"
     }
    }
  
  # return object of class Statmodel (from apply)
  .out  <- .StatModel(type = type, 
                      params = model,
                      varPosterior = as.numeric(NA),
                      dfPosterior = as.numeric(NA))
  return(.out)
    
  }, design = myDesign ) # end of apply here

  # Squeeze a set of sample variances together by computing empirical Bayes posterior means
  hlp <- limma::squeezeVar(var = sapply(models, getVar),
                           df = sapply(models, getDF))
  
  # put variance and degrees of freedom in appropriate slots
  for (i in 1:length(models)) {
    mydf <- hlp$df.prior + getDF(models[[i]])
    models[[i]]@varPosterior <- as.numeric(hlp$var.post[i])
    models[[i]]@dfPosterior <- as.numeric(mydf)
  }

  #return object of class StatModel
  return(models)
}
