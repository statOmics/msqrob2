#' @title The StatModel class for msqrob
#'
#' @description
#'
#' The `StatModel` class contains a statistical model as applied on a
#' feature.
#'
#' Models are created by the dedicated user-level functions
#' (`msqrob()`, `mqrobAggregate()`) or manually, using the
#' `StatModel()` constructor. In the former case, each quantitative
#' feature is assigned its statistical model and the models are stored
#' as a variable in a `DataFrame` object, as illustred in the example
#' below.
#'
#' @slot type `character(1)` defining type of the used model. Default
#'     is `"fitError"`, i.e. a error model. Other include `"lm"`,
#'     `"rlm"`, ...
#'
#' @slot params A `list()` containing information of the used model.
#'
#' @slot varPosterior `numeric()` of posterior variance.
#'
#' @slot dfPosterior `numeric()` of posterior degrees of freedom.
#'
#' @rdname StatModel
#'
#' @author Oliver M. Crook, Laurent Gatto, Lieven Clement
#'
#' @export
#'
#' @examples
#' ## A fully specified dummy model
#' myModel <- StatModel(
#'     type = "rlm",
#'     params = list(x = 3, y = 7, b = 4),
#'     varPosterior = c(0.1, 0.2, 0.3),
#'     dfPosterior = c(6, 7, 8)
#' )
#' myModel
#'
#'
#'
#' ## A collection of models stored as a variable in a DataFrame
#' mod1 <- StatModel(type = "rlm")
#' mod2 <- StatModel(type = "lm")
#' df <- DataFrame(x = 1:2)
#' df$mods <- c(mod1, mod2)
#' df
.StatModel <- setClass("StatModel",
    slots = c(
        type = "character",
        params = "list",
        varPosterior = "numeric",
        dfPosterior = "numeric"
    )
)

#' @rdname StatModel
#'
#' @param object `StatModel` object
#' @export
setMethod(
    "show", "StatModel",
    function(object) {
        cat("Object of class \"", class(object), "\"\n", sep = "")
        cat("The type is \"", object@type, "\"\n", sep = "")
        cat("There number of elements is", length(object@params), "\n")
        invisible(NULL)
    }
)


#' @rdname StatModel
#' @description Function for constructing a new `StatModel` object.
#'
#' @param type default set to fit-error, can be a "lm", "rlm" (robust lm with M estimation),
#'        "lmer" (when mixed models or ridge regression is adopted), "quasibinomial" (when peptide counts are fitted)
#'
#' @param params A list containing the parameters of the fitted model
#'
#' @param varPosterior Numeric, posterior variance, default is NA
#'
#' @param dfPosterior Numeric, posterior degrees of freedom,
#'     default is NA
#'
#' @examples
#' # TODO
#' @return A StatModel object
#'
#' @importFrom methods new
#'
#' @export
StatModel <- function(type = "fitError",
    params = list(),
    varPosterior = NA_real_,
    dfPosterior = NA_real_) {
    out <- new("StatModel")
    out@type <- type
    out@params <- params
    out@varPosterior <- varPosterior
    out@dfPosterior <- dfPosterior
    return(out)
}

setMethod("head", "StatModel", function(x, ...) x)

setMethod(
    "showAsCell", "StatModel",
    function(object) paste0("StatModel:", object@type)
)
