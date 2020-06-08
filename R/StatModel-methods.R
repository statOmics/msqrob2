setMethod("getContrast","StatModel",
          function(object, L){
              if (class(L) != "matrix") L <- as.matrix(L)
              coefs <- getCoef(object)
              out <- matrix(rep(NA, ncol(L)))
              rownames(out) <- colnames(L)
              hlp <- try(t(L) %*% coefs[rownames(L)], silent = TRUE)
              if (class(out)[1] != "try-error") out[]<-hlp
              return(out)
              })

setMethod("varContrast","StatModel",
          function(object, L){
              if (class(L) != "matrix") L <- as.matrix(L)
              out <- matrix(NA, ncol(L), ncol(L))
              rownames(out) <- colnames(out) <- colnames(L)
              vcovTmp <- getVcovUnscaled(object) * object@varPosterior
              hlp <- try(t(L) %*% vcovTmp[rownames(L), rownames(L)] %*% L,silent=TRUE)
              if (class(hlp)[1] != "try-error") out[]<-hlp
              return(out)
              })
