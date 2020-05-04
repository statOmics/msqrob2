#' Function to fit msqrob models
#'
#' @description Parameter estimation of msqrob models.
#'
#' @param y is a matrix with the quantified MS intensities, the features are in the rows and samples in the columns.
#' @param formula Model formula. The model is built based on the covariates in the data object.
#' @param data is a data frame of class DataFrame with information on the design. It has the same number of rows as the number of columns (samples) of y.
#' @param robust boolean to indicate if robust regression is performed to account for outliers. If 'FALSE' an OLS fit is performed.
#' @param maxitRob Maximum iterations in the IRWLS algorithm used in the M-estimation step of the robust regression.
#' @examples #TODO
#' @return A list of objects of the StatModel class
#' @rdname msqrob
#' @author Lieven Clement, Oliver M. Crook
#' @export
msqrobLm <- function(y,
                   formula,
                   data,
                   robust = TRUE,
                   maxitRob = 5)
{
  myDesign <- model.matrix(formula, data)
  models <- apply(y, 1, function(y, design)
  {
    # computatability check
    obs <- is.finite(y)
    type <- "fitError"
    model <- list(coefficients=NA,vcovUnscaled=NA,sigma=NA,df.residual=NA,w=NULL)

    if (sum(obs) > 0)  {
     # subset to finite observations, attention with R column switching
     X <- design[obs, , drop = FALSE]
     y <- y[obs]

     if(robust) {
     # use robust regression from MASS package, "M" estimation is used
     mod <- try(MASS::rlm(X,
                            y,
                            method = "M",
                            maxit = maxitRob),
                            silent = TRUE)
     if (class(mod)[1]!="try-error") type <- "rlm"
     } else {
     # if robust regression is not performed use standard linear fit
     mod <- try(lm.fit(X, y))
     if (class(mod)[1]!="try-error"&mod$rank==ncol(X)) type <- "lm"
     }

     if (type=="rlm") {
        w <- mod$w
        sigma <- sqrt(sum(mod$w*mod$resid^2)/(sum(mod$w)-mod$rank))
        df.residual<-sum(mod$w)-mod$rank
        }
     if (type=="lm") {
        w <- NULL
        sigma <- sqrt(sum(mod$residuals^2/mod$df.residual))
        df.residual<-mod$df.residual
        }
     if (type!="fitError") model<-list(coefficients=mod$coef,
        vcovUnscaled=.vcovUnscaled(mod),
        sigma=sigma,
        df.residual=df.residual,
        w=w)
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


#'@export
msqrobLmer <- function(y,formula,data,robust=TRUE,maxitRob=1,tol=1e-6,doQR=TRUE,featureGroups=NULL,lmerArgs = list(control = lmerControl(calc.derivs = FALSE)))
{
  require(lme4)

  if(length(formula)==3)
  formula<-formula[-2]

  fixed <- model.matrix(nobars(formula),data)

  df <- data
  df$fixed <- fixed

  if (ncol(fixed)<=2 & nobars(formula)[[2]]!=1)
  stop("Error: the mean model must have more than two parameters for ridge regression")

  if(is.null(findbars(formula)))
    form<-formula(y ~ (1|ridge)) else
    if(nobars(formula)[[2]]!=1)
      form<-formula(
        paste0("y ~ (1|ridge) + ",paste0("(",paste(findbars(formula),collapse=")+("),")"))
        ) else
      form <- formula

  if (is.null(featureGroups)) featureGroups<-rownames(y)
  y<-split.data.frame(y,featureGroups)

  models <- bplapply(y, function(y, form, data)
  {

    if(nrow(y)>1)
    {
      data<-data[rep(1:nrow(data),each=nrow(y)),]
      data$samples<-as.factor(rep(colnames(y),each=nrow(y)))
      data$features<-as.factor(rep(rownames(y),ncol(y)))
      form<-update.formula(form,~.+(1|samples)+(1|features))
    }
    dim(y)<-NULL
    data$y<-y
    data<-data[!is.na(data$y),]
    qrFixed<-qr(data$fixed)

    #for ridge regression of fixed effects
    if (doQR) {
      Q<-qr.Q(qrFixed)
    } else Q<-data$fixed

    model<-NULL

    if (nobars(formula)[[2]]!=1)
      ##Fooling lmer to adopt ridge regression using Fabian Scheipl's trick
      {if (qrFixed$rank==ncol(data$fixed)) try({
        data$ridge<-as.factor(rep(colnames(data$fixed)[-1],length=nrow(data)))
        parsedFormulaC <- lFormula(form,data=as.list(data))

        parsedFormulaC$reTrms$cnms$ridge <- ""
        ridgeId<-grep(names(parsedFormulaC$reTrms$Ztlist),pattern="ridge")
        colnames(Q)<-colnames(data$fixed)
        parsedFormulaC$reTrms$Ztlist[[ridgeId]]<-as(Matrix(t(Q[,-1])), class(parsedFormulaC$reTrms$Ztlist[[ridgeId]]))
        parsedFormulaC$reTrms$Zt<-do.call(rbind,parsedFormulaC$reTrms$Ztlist)

        devianceFunctionC <- do.call(mkLmerDevfun, parsedFormulaC)
        optimizerOutputC <- optimizeLmer(devianceFunctionC)
        model <- mkMerMod(
                         rho = environment(devianceFunctionC),
                         opt = optimizerOutputC,
                         reTrms = parsedFormulaC$reTrms,
                        fr = parsedFormulaC$fr)
        },silent=TRUE)} else
        ##For advanced users who opted to perform ridge regression by specifying all effects as random effects
        try(
        model <-lmer(form,as.list(data))
        ,silent=TRUE)

    if (is.null(model)) {
      type <- "fitError"
      model <- list(coefficients=NA,vcovUnscaled=NA,sigma=NA,df.residual=NA)
    } else {
      type <- "lmer"
      sseOld <- model@devcomp$cmp['pwrss']
      while (maxitRob > 0){
        maxitRob <- maxitRob-1
        res <- resid(model)
        model@frame$`(weights)` <- MASS::psi.huber(res/(mad(res,0)))
        model <- refit(model)
        sse <- model@devcomp$cmp['pwrss']
        if(abs(sseOld-sse)/sseOld <= tol) break
        sseOld <- sse
      }

      sigma<-sigma(model)
      betas <-.getBetaB(model)
      vcovUnscaled<-as.matrix(.getVcovBetaBUnscaled(model))
      if (nobars(formula)[[2]]!=1)
      {
        if (doQR)
        {
          ids<-c(1,grep("ridge",names(betas)))
          Rinv <- diag(betas)
          coefNames<-names(betas)
          Rinv[ids,ids]<-solve(qr.R(qrFixed))
          Rinv[1, 1] <- 1
          betas <- c(Rinv%*%betas)
          names(betas)<-coefNames
          vcovUnscaled<-t(Rinv)%*%vcovUnscaled%*%Rinv
          rownames(vcovUnscaled) <- colnames(vcovUnscaled) <- names(betas)
        }
      }
      df.residual<-.getDfLmer(model)
      if (df.residual<2L)
        model<-list(coefficients=NA,vcovUnscaled=NA,sigma=NA,df.residual=NA,w=NA) else
        model<-list(coefficients=betas,vcovUnscaled=vcovUnscaled,sigma=sigma,df.residual=df.residual,w=model@frame$`(weights)`)
    }
  return(StatModel(type=type,params=model,varPosterior=as.numeric(NA),dfPosterior=as.numeric(NA)))
  },form=update.formula(form,y~.),data=df)

  hlp<-limma::squeezeVar(var=sapply(models,getVar),df=sapply(models,getDF))
  for (i in 1:length(models))
  {
    models[[i]]@varPosterior<-as.numeric(hlp$var.post[i])
    models[[i]]@dfPosterior<-as.numeric(hlp$df.prior+getDF(models[[i]]))
  }
  return(models)
}

## Calculate unscaled covariance matrix for lm or rlm fit
.vcovUnscaled <-function(model)
{
p1 <- 1L:model$rank
p <- length(model$coefficients)
out<-matrix(NA,p,p)
out[p1,p1]<-chol2inv( model$qr$qr[p1, p1, drop = FALSE])
colnames(out)<-rownames(out)<-names(model$coefficients)
return(out)
}

#' @import lme4
.getVcovBetaBUnscaled <- function(model){
  X <- getME(model,"X")
  Z <- getME(model,"Z")
  vcovInv <- Matrix::crossprod(cbind2(X,Z))
  Ginv <- Matrix::solve(Matrix::tcrossprod(getME(model,"Lambda")) +
                        Matrix::Diagonal(ncol(Z),1e-18))
  i <- -seq_len(ncol(X))
  vcovInv[i,i] <- vcovInv[i,i]+Ginv
  vcovInv<-Matrix::solve(vcovInv)
  ranefLevels<-imap(model@flist,~{paste0(.y,levels(.x))})
  zNames<-unlist(lapply(1:length(model@cnms),function(x,cnms,levels) c(outer(cnms[[x]],levels[[names(cnms)[x]]],paste0)),cnms=model@cnms,levels=ranefLevels))


  rownames(vcovInv)<-colnames(vcovInv)<-c(colnames(X), zNames)
  return(vcovInv)
}

#' @import purrr
#' @import lme4
.getBetaB <- function(model) {
  betaB <- c(as.vector(getME(model,"beta")),as.vector(getME(model,"b")))
  ranefLevels<-imap(model@flist,~{paste0(.y,levels(.x))})
  zNames<-unlist(lapply(1:length(model@cnms),function(x,cnms,levels) c(outer(cnms[[x]],levels[[names(cnms)[x]]],paste0)),cnms=model@cnms,levels=ranefLevels))
  names(betaB) <- c(colnames(model@pp$X), zNames)
  betaB
}

#' @import lme4
.getDfLmer <- function(object){
  w <- object@frame$"(weights)"
  if (is.null(w)) w <- 1
  sigma <- sigma(object)
  sum((resid(object)* sqrt(w))^2)/sigma^2
}
