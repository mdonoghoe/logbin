#  Based on File src/library/stats/R/predict.glm.R
#       and File src/library/stats/R/lm.R
#  Part of the R package, https://www.R-project.org
# 
#  Modified by Mark W. Donoghoe
#     28/08/2015 - to work with logbin objects
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

predict.logbin <- function(object, newdata = NULL, type = c("link", "response",
    "terms"), terms = NULL, na.action = na.pass, checkminmax = TRUE,...)
{
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    
    tt <- terms(object)
    if (missing(newdata) || is.null(newdata)) {
        mm <- X <- model.matrix(object)
        mmDone <- TRUE
        offset <- object$offset
    } else {
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action, xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)
        if (checkminmax)
            .checkXminmax(m, xminmax = object$xminmax)
        X <- model.matrix(Terms, m)
        offset <- rep(0, nrow(X))
        if (!is.null(off.num <- attr(tt, "offset")))
            for (i in off.num) offset <- offset + eval(attr(tt, "variables")[[i+1]], newdata)
        if (!is.null(object$call$offset))
            offset <- offset + eval(object$call$offset, newdata)
        mmDone <- FALSE
    }
    n <- length(object$residuals)
    p <- object$rank
    p1 <- seq_len(p)
    beta <- object$coefficients
    predictor <- drop(X %*% beta)
    if(!is.null(offset))
        predictor <- predictor + offset
    if (type == "terms") {
        if (!mmDone) {
            mm <- model.matrix(object)
            mmDone <- TRUE
        }
        aa <- attr(mm, "assign")
        ll <- attr(tt, "term.labels")
        hasintercept <- attr(tt, "intercept") > 0L
        if (hasintercept)
            ll <- c("(Intercept)", ll)
        aaa <- factor(aa, labels = ll)
        asgn <- split(order(aa), aaa)
        if (hasintercept) {
            asgn$"(Intercept)" <- NULL
            if(!mmDone) {
                mm <- model.matrix(object)
                mmDone <- TRUE
            }
            avx <- colMeans(mm)
            termsconst <- sum(avx * beta)
        }
        nterms <- length(asgn)
        if (nterms > 0) {
            predictor <- matrix(ncol = nterms, nrow = NROW(X))
            dimnames(predictor) <- list(rownames(X), names(asgn))
            if (hasintercept)
                X <- sweep(X, 2L, avx, check.margin = FALSE)
            unpiv <- rep.int(0L, NCOL(X))
            unpiv <- p1
            for (i in seq.int(1L, nterms, length.out = nterms)) {
                iipiv <- asgn[[i]]
                ii <- unpiv[iipiv]
                iipiv[ii == 0L] <- 0L
                predictor[, i] <- if (any(iipiv > 0L))
                    X[, iipiv, drop = FALSE] %*% beta[iipiv]
                else 0
            }
            if (!is.null(terms))
                predictor <- predictor[, terms, drop = FALSE]
        } else {
            predictor <- matrix(0, n, 0L)
        }
        attr(predictor, "constant") <- if (hasintercept)
            termsconst
        else 0
    }
    if (missing(newdata) && !is.null(na.act))
        predictor <- napredict(na.act, predictor)
   
    switch(type, response = {predictor <- family(object)$linkinv(predictor)},
            link = , terms =)
    predictor
}