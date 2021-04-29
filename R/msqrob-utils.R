#' Smallest unique protein groups
#'
#' @description For a given vector of protein group names, outputs the names of those protein groups for which none of its member proteins is present in a smaller protein group.
#' @param proteins A vector of characters or factors containing single proteins and/or protein groups (i.e. proteins separated by a separator symbol).
#' @param split The character string that is used to separate the indivudual protein names in each protein group.
#' @return A character vector containing the names of the protein groups for which none of its proteins is present in a smaller protein group.
#' @examples
#' data(pe)
#' smallestUniqueGroups(rowData(pe[["peptide"]])$Proteins)
#' @export
smallestUniqueGroups <- function(proteins,
    split = ";") {
    b <- strsplit(x = as.character(unique(proteins)), split = split, fixed = TRUE)

    included <- vector()

    j <- 1
    while (length(b) != 0) {
        included <- c(
            included,
            vapply(
                b[vapply(b, length, integer(1)) == j],
                function(x) paste(x, collapse = split),
                character(1)
            )
        )
        a <- unlist(b[vapply(b, length, integer(1)) == j])
        b <- b[vapply(b, length, integer(1)) > j]

        if (length(b) != 0) {
            sel <- vector()
            for (i in seq_len(length(b))) {
                sel[i] <- !any(b[[i]] %in% a)
            }
            b <- b[sel]
            j <- j + 1
        }
    }

    included <- unlist(included)
    return(included)
}

#' Make contrast matrix
#'
#' @description  Construct the contrast matrix corresponding to specified contrasts
#'               of a set of parameters.
#'
#' @param contrasts character vector specifying contrasts, i.e. the linear combination of the modelparameters that equals to zero.
#'
#' @param parameterNames character vector specifying the model parameters that are involved in the contrasts, e.g if we model data of
#'        three conditions using a factor condition with three levels a, b and c then our model will have 3 mean parameters named (Intercept),
#'        conditionb and conditionc. Hence the log2 fold change between  b and a is conditionb. Under the null hypothesis the log2 fold change
#'        equals 0. Which is to be encoded as "conditionb=0". If we would like to test for log2 fold change between condition c and b we assess if
#'        the log2 fold change conditionc-conditionb equals 0, encoded as "condtionb-conditionc=0".
#'
#' @examples
#' makeContrast(c("conditionb = 0"),
#'     parameterNames = c(
#'         "(Intercept)",
#'         "conditionb",
#'         "conditionc"
#'     )
#' )
#' makeContrast(c("conditionc=0"),
#'     parameterNames = c("conditionc")
#' )
#' makeContrast(c(
#'     "conditionb=0",
#'     "conditionc=0",
#'     "conditionc-conditionb=0"
#' ),
#' parameterNames = c(
#'     "conditionb",
#'     "conditionc"
#' )
#' )
#' @return A numeric contrast matrix with rownames that equal the model parameters that are involved in the contrasts
#'
#' @rdname makeContrast
#'
#' @export


makeContrast <- function(contrasts, parameterNames) {
    return(t(.chrlinfct2matrix(contrasts, parameterNames)$K))
}

##### None exported functions from multcomp package is included here to
##### During R and Bioc checks


.chrlinfct2matrix <- function(ex, var) {
    if (!is.character(ex)) {
        stop("msqrob2:::.chrlinfct2matrix: argument ", sQuote(ex),
            " is not of type character",
            call. = FALSE
        )
    }
    if (!is.character(var)) {
        stop("msqrob2:::.chrlinfct2matrix: argument ", sQuote(var),
            " is not of type character",
            call. = FALSE
        )
    }
    K <- matrix(0, nrow = length(ex), ncol = length(var))
    colnames(K) <- var
    rownames(K) <- seq_along(ex)
    m <- rep(0, length(ex))
    for (i in seq_len(length(ex))) {
        expr <- parse(text = ex[i])
        if (length(expr[[1]]) != 3) {
            stop("msqrob2:::.chrlinfct2matrix: argument ", sQuote(ex[i]),
                " cannot be interpreted as expression",
                call. = FALSE
            )
        }
        tmp <- .expression2coef(expr, vars = var)
        if ("(Intercept)" %in% var) {
            tmp$names[tmp$names == "Intercept"] <- "(Intercept)"
        }
        if (!all(tmp$names %in% var)) {
            stop("msqrob2:::.chrlinfct2matrix: variable(s) ",
                paste(sQuote(tmp$names[!tmp$names %in% var]),
                    collapse = ", "
                ), " not found",
                call. = FALSE
            )
        }
        for (n in tmp$names) {
            K[i, var == n] <- tmp$coef[tmp$names ==
                n]
        }
        m[i] <- tmp$m
        if (i == 1) {
            alternative <- tmp$alternative
        }
        else {
            if (tmp$alternative != alternative) {
                stop("msqrob2:::.chrlinfct2matrix: mix of alternatives currently not implemented",
                    call. = FALSE
                )
            }
        }
        rownames(K)[i] <- paste0(tmp$lhs, collapse = "")
    }
    list(K = K, m = m, alternative = alternative)
}

#' @importFrom codetools walkCode makeCodeWalker
#' @importFrom methods is

.expression2coef <- function(ex, vars, debug = FALSE) {
    m.rhs <- .rhs(ex)
    m.lhs <- .lhs(ex)
    symcoef <- new.env(parent = emptyenv())
    get_coef_attr <- function(x) {
        if (is.symbol(x)) {
            get0(as.character(x), envir = symcoef, ifnotfound = NULL)
        } else {
            attr(x, "coef")
        }
    }
    set_coef_attr <- function(x, val) {
        if (is.symbol(x)) {
            if (is.null(x)) {
                rm(as.character(x), envir = symcoef)
            } else {
                assign(as.character(x), val, envir = symcoef)
            }
        }
        else {
            attr(x, "coef") <- val
        }
        x
    }
    set_coef_attr(m.lhs, 1)
    if (debug) {
        message(".expression2coef", ": lhs is ", sQuote(paste0(deparse(m.lhs),
            collapse = ""
        )))
        message(".expression2coef", ": rhs is ", sQuote(paste0(deparse(m.rhs),
            collapse = ""
        )))
    }
    effects <- walkCode(m.lhs, makeCodeWalker(handler = function(v,
    w) {
        if (debug) {
            w$trace("handler", v, w)
        }
        switch(v,
            `-` = w$sub,
            `+` = w$add,
            `*` = w$mul,
            `/` = w$div,
            `(` = w$exp,
            `:` = w$ita,
            w$eval
        )
    }, is.effect <- function(x) {
        as.character(x) %in% vars
    }, eval = function(v, w) {
        if (debug) {
            w$trace("eval", v, w)
        }
        parms <- c()
        for (e in as.list(v)[-1]) {
            coef <- 1
            if (debug) {
                message(
                    "eval", ": walking ", sQuote(e), " with coef = ",
                    coef
                )
            }
            parms <- c(parms, p <- walkCode(
                w$setCoef(e, coef),
                w
            ))
            if (is.effect(p)) {
                w$fatal(
                    "eval", "within ", sQuote(deparse(v)),
                    ", the term ", sQuote(p), " ", "must not denote an effect. Apart from that, ",
                    "the term must evaluate to a real valued constant"
                )
            }
        }
        cparms <- c()
        for (e in parms) {
            cparms <- c(cparms, parse(text = paste(e, "*", w$getCoef(e))))
        }
        if (debug) {
            dumped <- lapply(cparms, function(x, w) {
                paste(
                    x,
                    "with coef =", w$getCoef(x)
                )
            }, w)
            message("eval", ": cparms = ", w$enum(dumped))
        }
        res <- try(do.call(as.character(v[[1]]), as.list(cparms)),
            silent = TRUE
        )
        if (is(res, "try-error")) {
            w$fatal(
                "eval", "the evaluation of the expression ",
                sQuote(deparse(v)), " ", "failed with ", dQuote(attr(
                    res,
                    "condition"
                )$message)
            )
        }
        if (length(res) != 1 || !is.numeric(res) || !is.finite(res)) {
            w$fatal(
                "eval", "the expression ", sQuote(deparse(v)),
                " ", "did not evaluate to a real valued constant. ",
                "Result is ", sQuote(res)
            )
        }
        res <- w$setCoef(res * w$getCoef(v), 1)
        if (debug) {
            dumped <- lapply(res, function(x, w) {
                paste(
                    x, "with coef =",
                    w$getCoef(x)
                )
            }, w)
            message("eval", ": res = ", w$enum(dumped))
        }
        res
    }, call = function(v, w) {
        if (debug) {
            w$trace("call", v, w)
        }
        w$fatal(
            "call", "there is probably a syntax error within subexpression",
            sQuote(deparse(v))
        )
    }, sub = function(v, w) {
        if (debug) {
            w$trace("sub", v, w)
        }
        uminus <- length(as.list(v)) == 2
        minuend <- as.list(v)[2]
        subtrahend <- as.list(v)[3]
        if (debug) {
            message(
                "sub", ": minuend is ", sQuote(minuend),
                ", coef = ", w$getCoef(minuend)
            )
            message(
                "sub", ": subtrahend is ", sQuote(subtrahend),
                ", coef = ", w$getCoef(subtrahend)
            )
            message("sub", ": uminus is ", uminus)
        }
        exp.coef <- ifelse(uminus, -w$getCoef(v), w$getCoef(v))
        res <- c()
        for (e in minuend) {
            if (debug) {
                message(
                    "sub", ": walking minuend ", sQuote(e),
                    ", coef = ", exp.coef
                )
            }
            res <- c(res, walkCode(w$setCoef(e, exp.coef), w))
        }
        if (!uminus) {
            for (e in subtrahend) {
                if (debug) {
                    message(
                        "sub", ": walking subtrahend ", sQuote(e),
                        ", coef = ", -exp.coef
                    )
                }
                res <- c(res, walkCode(
                    w$setCoef(e, -exp.coef),
                    w
                ))
            }
        }
        sum <- 0
        symbols <- c()
        for (e in res) {
            if (is.numeric(e)) {
                sum <- sum + e
            } else {
                symbols <- c(symbols, e)
            }
        }
        if (length(dups <- symbols[duplicated(symbols)])) {
            w$fatal(
                "sub", "multiple occurence of ", w$enum(dups),
                " ", "found within expression ", sQuote(deparse(v))
            )
        }
        if (length(symbols) == 0) {
            return(w$setCoef(sum, 1))
        }
        if (sum) {
            w$fatal(
                "sub", "forming a difference between a constant and ",
                "an effect as in ", sQuote(deparse(v)), " ",
                "is not supported"
            )
        }
        symbols
    }, ita = function(v, w) {
        if (debug) {
            w$trace("ita", v, w)
        }
        tmp <- deparse(v)
        prefix <- gsub("^([+-]*)(.*)", "\\1", tmp)
        name <- gsub("^([+-]*)(.*)", "\\2", tmp)
        sign <- 1
        if (prefix != "") {
            for (x in base::unlist(base::strsplit(prefix, ""))) {
                sign <- sign * switch(x,
                    `-` = -1,
                    `+` = 1,
                    w$fatal(
                        "ita",
                        "strange character ", sQuote(x), " seen in ",
                        sQuote(deparse(v))
                    )
                )
            }
        }
        res <- w$setCoef(as.name(name), sign * w$getCoef(v))
        if (debug) {
            dumped <- lapply(res, function(x, w) {
                paste(
                    x, "with coef =",
                    w$getCoef(x)
                )
            }, w)
            message("ita", ": res = ", w$enum(dumped))
        }
        res
    }, exp = function(v, w) {
        if (debug) {
            w$trace("exp", v, w)
        }
        res <- c()
        for (e in as.list(v)[-1]) {
            res <- c(res, walkCode(w$setCoef(e, 1), w))
        }
        symbols <- c()
        for (e in res) {
            symbols <- c(symbols, w$setCoef(e, w$getCoef(e) *
                w$getCoef(v)))
        }
        if (debug) {
            dumped <- lapply(symbols, function(x, w) {
                paste(
                    x,
                    "with coef =", w$getCoef(x)
                )
            }, w)
            message("exp", ": res = ", w$enum(dumped))
        }
        symbols
    }, add = function(v, w) {
        if (debug) {
            w$trace("add", v, w)
        }
        res <- c()
        for (e in as.list(v)[-1]) {
            res <- c(res, walkCode(w$setCoef(e, 1), w))
        }
        symbols <- c()
        sum <- 0
        for (e in res) {
            if (is.numeric(e)) {
                sum <- sum + e
            } else {
                symbols <- c(symbols, e)
            }
        }
        if (length(symbols) == 0) {
            return(w$setCoef(sum * w$getCoef(v), 1))
        }
        if (length(dups <- symbols[duplicated(symbols)]) != 0) {
            w$fatal(
                "add", "multiple occurence of ", w$enum(dups),
                " ", "within subexpression ", sQuote(deparse(v))
            )
        }
        if (sum) {
            w$fatal(
                "add", "adding up a constant and an effect ",
                "as in ", sQuote(deparse(v)), " is not supported"
            )
        }
        res <- c()
        for (e in symbols) {
            res <- c(res, w$setCoef(e, w$getCoef(e) * w$getCoef(v)))
        }
        res
    }, mul = function(v, w) {
        if (debug) {
            w$trace("mul", v, w)
        }
        res <- c()
        for (e in as.list(v)[-1]) {
            res <- c(res, walkCode(w$setCoef(e, 1), w))
        }
        product <- 1
        symbols <- c()
        for (r in res) {
            if (is.numeric(r)) {
                product <- product * r
            } else {
                symbols <- c(symbols, r)
            }
        }
        if (product == 0 && length(symbols)) {
            w$fatal(
                "mul", "The constant part of the expression ",
                sQuote(deparse(v)), " ", "evaluates to zero. This would zero out the effect(s) ",
                sQuote(symbols)
            )
        }
        product <- product * w$getCoef(v)
        if (length(symbols) == 0) {
            return(w$setCoef(product, 1))
        }
        if (length(symbols) > 1 && all(unlist(lapply(res, is.symbol)))) {
            w$fatal(
                "mul", "the multiplication of effects ",
                w$enum(symbols), " ", "as in ", sQuote(deparse(v)),
                " is not supported"
            )
        }
        res <- c()
        for (s in symbols) {
            res <- c(res, w$setCoef(s, w$getCoef(s) * product))
        }
        if (debug) {
            dumped <- lapply(res, function(x, w) {
                paste(
                    x, "with coef =",
                    w$getCoef(x)
                )
            }, w)
            message("mul", ": res = ", w$enum(dumped))
        }
        res
    }, div = function(v, w) {
        if (debug) {
            w$trace("div", v, w)
        }
        if ((lv <- length(v)) != 3) {
            w$fatal(
                "div", "internal error: length of language object ",
                sQuote(v), " ", "is not 3, but ", lv, ". Please file a bug report"
            )
        }
        dividend <- c()
        for (e in as.list(v)[2]) {
            dividend <- c(dividend, walkCode(
                w$setCoef(e, 1),
                w
            ))
        }
        divisor <- c()
        for (e in as.list(v)[3]) {
            divisor <- c(divisor, walkCode(w$setCoef(e, 1), w))
        }
        if (length(divisor) != 1) {
            w$fatal(
                "div", "can't divide by ", sQuote(divisor),
                " in ", sQuote(deparse(v))
            )
        }
        if (any(unlist(lapply(divisor, is.effect)))) {
            w$fatal(
                "div", "cant't divide by effect ", sQuote(divisor),
                " in ", sQuote(deparse(v))
            )
        }
        if (any(unlist(lapply(divisor, is.symbol)))) {
            w$fatal(
                "div", "cant't divide by symbol ", sQuote(divisor),
                " in ", sQuote(deparse(v))
            )
        }
        divisor <- as.numeric(divisor)
        if (!is.finite(divisor) || divisor == 0) {
            w$fatal(
                "div", "can't divide by ", sQuote(divisor),
                " in ", sQuote(deparse(v))
            )
        }
        res <- c()
        for (s in dividend) {
            if (is.numeric(s)) {
                res <- c(res, s * w$getCoef(v) / divisor)
            }
            else {
                res <- c(res, w$setCoef(s, w$getCoef(s) * w$getCoef(v) / divisor))
            }
        }
        if (debug) {
            message("div", ": dividend = ", w$enum(dividend))
            message("div", ": divisor = ", w$enum(divisor))
            message("div", ": res = ", w$enum(res))
        }
        res
    }, leaf = function(e, w) {
        if (debug) {
            w$trace("leaf", e, w)
        }
        if (is.numeric(e)) {
            return(w$setCoef(e * w$getCoef(e), 1))
        }
        e
    }, getCoef = function(e) {
        a <- get_coef_attr(e)
        ifelse(is.null(a), 1, a)
    }, setCoef = function(e, coef) {
        set_coef_attr(e, coef)
    }, enum = function(x) {
        paste0("'", x, "'", collapse = ", ")
    }, fatal = function(name, ...) {
        stop(
            "msqrob2:::.expression2coef::walkCode::",
            name, ": ", ...,
            call. = FALSE
        )
    }, trace = function(fn, v, w) {
        message(
            fn, ": v = ", sQuote(v), ", mode = ", mode(v),
            ", typeof = ", typeof(v), ", length = ", length(v),
            ", coef   = ", w$getCoef(v)
        )
    }))
    if (any(idx <- is.numeric(effects))) {
        stop("msqrob2:::.expression2coef: The lhs expression ",
            sQuote(deparse(m.lhs)), " ", "contains a numeric offset term evaluating to ",
            paste0(effects[idx], collapse = ", "), ". ", "This is either an internal flaw or a misspecification from your part. ",
            "If so, please pull these offsets to the right-hand side of the equation",
            call. = FALSE
        )
    }
    effect.names <- c()
    effect.coefs <- c()
    for (effect in c(effects)) {
        effect.names <- c(effect.names, as.character(effect))
        effect.coefs <- c(effect.coefs, get_coef_attr(effect))
    }
    list(
        coef = effect.coefs, names = effect.names, m = m.rhs,
        alternative = .side(ex), lhs = deparse(m.lhs, width.cutoff = 500)
    )
}


.rhs <- function(ex) {
    if (length(ex) != 1) {
        stop("msqrob2:::.rhs: expression is not of length 1",
            call. = FALSE
        )
    }
    if (length(ex[[1]][[3]]) == 2) {
        return(-ex[[1]][[3]][[2]])
    }
    rhs <- ex[[1]][[3]]
    if (!.is_num(rhs) || length(rhs) > 1) {
        stop("msqrob2:::.rhs: right hand side of expression ",
            sQuote(ex), " is not a scalar numeric",
            call. = FALSE
        )
    }
    return(rhs)
}

.lhs <- function(ex) {
    if (length(ex) != 1) {
        stop("msqrob2:::.lhs: expression is not of length 1",
            call. = FALSE
        )
    }
    if (length(ex[[1]]) != 3) {
        stop("msqrob2:::.lhs: expression ", sQuote(ex), " does not contain a left and right hand side",
            call. = FALSE
        )
    }
    return(ex[[1]][[2]])
}

.is_num <- function(x) {
    if (length(x) == 1) {
          return(is.numeric(x))
      }
    if (length(x) == 2) {
          return(is.name(x[[1]]) && is.numeric(x[[2]]))
      }
    return(FALSE)
}

.side <- function(ex) {
    side <- .as.char(ex[[1]][[1]])
    if (!(side %in% c("<=", ">=", "==", "="))) {
        stop("msqrob2:::side: does not contain ", sQuote("<=, >=, =="),
            call. = FALSE
        )
    }
    alternative <- switch(side,
        `<=` = "greater",
        `>=` = "less",
        `==` = "two.sided",
        `=` = "two.sided"
    )
    return(alternative)
}

.as.char <- function(ex) {
    if (length(ex) == 1) {
          return(as.character(ex))
      }
    if (length(ex) == 3 && ex[[1]] == ":") {
          return(paste(.as.char(ex[[2]]), ":", .as.char(ex[[3]]),
              sep = ""
          ))
      }
    stop(
        "msqrob2:::.as.char: Failed to convert expression ",
        ex, " to character"
    )
}
