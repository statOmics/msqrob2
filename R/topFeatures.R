#' Toplist of DE proteins, peptides or features
#'
#' @description Summary table of the differentially expressed Features
#'
#' @param models A list with elements of the class `StatModel` that are
#'        estimated using the \code{\link{msqrob}} function
#' @param contrast `numeric` (matrix)vector specifying one contrast of
#'        the linear model coefficients to be tested equal to zero.
#'        The (row)names of the vector should be equal to the names of
#'        parameters of the model.
#' @param adjust.method `character` specifying the method to adjust
#'        the p-values for multiple testing.
#'        Options, in increasing conservatism, include ‘"none"’,
#'        ‘"BH"’, ‘"BY"’ and ‘"holm"’.  See ‘p.adjust’ for the complete
#'        list of options. Default is "BH" the Benjamini-Hochberg method
#'        to control the *False Discovery Rate* (FDR).
#' @param ddf.method `character` specifying the method to calculate the
#'       *denominator degrees of freedom* (DoF) for the t-test.
#'       Options are:
#'         - "residual": (default) use the posterior residual DoF (\code{getDfPosterior}),
#'           which was the only DDoF method for *msqrob2* before v2.XX.
#'           It can significantly overestimate the denominator DoF,
#'           resulting in overly significant p-values.
#'         - "ML1": use the `dof_ml1` function from the *parameters* package for
#'           heuristic approximation of the DDoF, which is much faster than
#'           the Kenward-Roger method, while still providing an adequate approximation.
#'         - "KenwardRoger": use the `dof_kenward` function from the *parameters* package for
#'           *Kenward-Roger approximation* of the DDoF. This is the most accurate, but
#'           also most computationally intensive method.
#'         - "Satterthwaite": use the `dof_satterthwaite` function from the *parameters* package for
#'           *Satterthwaite approximation* of the DDoF. This is a faster, but less accurate
#'           alternative to *Kenward-Roger* method, which may not work well for small sample sizes.
#' @param sort `boolean(1)` to indicate if the features have to be sorted according
#'        to statistical significance.
#' @param alpha `numeric` specifying the cutoff value for adjusted p-values.
#'        Only features with lower p-values are listed.
#'
#' @examples
#' data(pe)
#'
#' # Aggregate peptide intensities in protein expression values
#' pe <- aggregateFeatures(pe, i = "peptide", fcol = "Proteins", name = "protein")
#'
#' # Fit msqrob model
#' pe <- msqrob(pe, i = "protein", formula = ~condition)
#'
#' # Define contrast
#' getCoef(rowData(pe[["protein"]])$msqrobModels[[1]])
#'
#' # Assess log2 fold change between condition c and condition b:
#' L <- makeContrast("conditionc - conditionb=0", c("conditionb", "conditionc"))
#' topDeProteins <- topFeatures(rowData(pe[["protein"]])$msqrobModels, L)
#' @return A dataframe with log2 fold changes (logFC), standard errors (se),
#'         degrees of freedom of the test (df), t-test statistic (t),
#'         p-values (pval) and adjusted pvalues (adjPval) using the specified
#'         adjust.method in the p.adjust function of the stats package.
#' @importFrom stats pt p.adjust na.exclude
#' @rdname topTable
#' @author Lieven Clement
#' @export
topFeatures <- function(
    models,
    contrast,
    adjust.method = "BH",
    ddf.method = c("residual", "ML1", "KenwardRoger", "Satterthwaite"),
    sort = TRUE,
    alpha = 1
) {
    if (!is(contrast, "matrix")) {
        contrast <- as.matrix(contrast)
    }
    if (ncol(contrast) > 1) {
        stop("The 'contrast' argument has more than one column, only one contrast is allowed")
    }
    # remove unused coefficients
    contrast <- contrast[rowSums(contrast) != 0, , drop = FALSE]

    logFC <- vapply(models,
        getContrast,
        numeric(1),
        L = contrast
    )
    se <- sqrt(vapply(models,
        varContrast,
        numeric(1),
        L = contrast
    ))
    ddf.method <- match.arg(ddf.method)
    if (ddf.method == "residual") {
        df <- vapply(models, getDfPosterior, numeric(1))
    } else if (ddf.method %in% c("KenwardRoger", "ML1", "Satterthwaite")) {
        if (!requireNamespace("parameters", quietly = TRUE)) {
            stop("parameters package is required to calculate the denominator DoF using ", ddf.method, " method.")
        }
        ddf_func <- if (ddf.method == "KenwardRoger") {
            parameters::dof_kenward
        } else if (ddf.method == "ML1") {
            parameters::dof_ml1
        } else if (ddf.method == "Satterthwaite") {
            parameters::dof_satterthwaite
        } else {
            stop("Unsupported ddf.method=", ddf.method)
        }
        df <- bplapply(models, function(model) {
            if ("model" %in% names(model@params) && is(model@params$model, "lmerMod")) {
                lmm <- model@params$model
                tryCatch(min(ddf_func(lmm)[rownames(contrast)]),
                            error = function(e) NA_real_)
            } else {
                NA_real_
            }
        })
        df <- unlist(df, recursive = FALSE, use.names = FALSE)
    } else {
        stop("Unsupported ddf.method=", ddf.method)
    }
    t <- logFC / se
    pval <- pt(-abs(t), df) * 2
    adjPval <- p.adjust(pval, method = adjust.method)
    out <- data.frame(logFC, se, df, t, pval, adjPval)
    rownames(out) <- names(models)
    if (alpha < 1) {
        signif <- adjPval < alpha
        out <- na.exclude(out[signif, ])
    }
    if (sort) {
        out <- out[order(out$pval), ]
    }
    return(out)
}
