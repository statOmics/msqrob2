##' classes for msqrob
##'
##' @title The StatModel class for msqrob
##' @slot type The type of the used model
##' @slot params A list containing information of the used model
##' @slot varPosterior A vector of posterior variance
##' @slot dfPosterior A vecotr of posterior degrees of freedom
##' @rd
##' @rdname msqrob-framework
##' @author Oliver M. Crook, Lieven Clement
##' @examples
##' myModel <- new("StatModel", type = "rlm", params = list(x = 3, y = 7, b = 4), varPosterior = c(0.1, 0.2, 0.3), dfPosterior = c(6, 7, 8))
##' myModel@params
##' mod1 <- new("StatModel", type = "rlm")
##' mod2 <- new("StatModel", type = "lm")
##' df <- DataFrame(x = 1:2)
##' df$mods <- c(mod1, mod2)
##' df$mods
#' @export
setClass("StatModel",
         slots = c(type = "character",
                   params = "list",
                   varPosterior = "numeric",
                   dfPosterior = "numeric"))

##' @title show method for object of class StatModel
##'
##' @rdname msqrob-framework
setMethod("show", "StatModel",
          function(object) {
            cat("Object of class \"", class(object), "\"\n", sep = "")
            cat("The type is \"", object@type, "\"\n", sep = "")
            cat("There number of elements is", length(object@params),"\n")
            invisible(NULL)
          })
