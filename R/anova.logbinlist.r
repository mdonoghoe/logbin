#  Based on File src/library/stats/R/glm.R
#  Part of the R package, https://www.R-project.org
#
#  Modified by Mark W. Donoghoe
#     28/08/2015 - to work with logbin objects
#
#  Copyright (C) 1995-2020 The R Core Team
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

anova.logbinlist <- function(object, ..., test = NULL) {
  doscore <- !is.null(test) && test == "Rao"
  responses <- as.character(lapply(object, function(x) {deparse(formula(x)[[2L]])}))
  sameresp <- responses == responses[1L]
  if (!all(sameresp)) {
    object <- object[sameresp]
    warning(gettextf("models with response %s removed because response differs from model 1", 
            sQuote(deparse(responses[!sameresp]))), domain = NA)
  }
  ns <- sapply(object, function(x) length(x$residuals))
  if (any(ns != ns[1L]))
    stop("models were not all fitted to the same size of dataset")
  boundary <- sapply(object, "[[", "boundary")
  if(any(boundary)) {
    object <- object[!boundary]
    warning(paste0("model(s) ", paste(which(boundary), collapse = ", ")," removed because MLE is on the boundary of the parameter space"))
  }
  nmodels <- length(object)
  if (nmodels == 0)
    stop("anova.logbinlist: no models remaining")
  if(nmodels == 1)
    stop("anova.logbin does not support anova for a single model. Fit nested models manually and input to anova.logbin")
  resdf <- as.numeric(lapply(object, function(x) x$df.residual))
  resdev <- as.numeric(lapply(object, function(x) x$deviance))
  if (doscore) {
    score <- numeric(nmodels)
    score[1] <- NA
    df <- -diff(resdf)
    for(i in seq_len(nmodels - 1)) {
      m1 <- if(df[i] > 0)
        object[[i]]
      else object[[i + 1]]
      m2 <- if(df[i] > 0)
        object[[i + 1]]
      else object[[i]]
      r <- m1$residuals
      w <- m1$prior.weights * m1$family$mu.eta(m1$linear.predictors)^2/m1$family$variance(m1$fitted.values)
      zz <- glm.fit(x = m2$x, y = r, weights = w)
      score[i + 1] <- zz$null.deviance - zz$deviance
      if (df[i] < 0) score[i + 1] <- -score[i + 1]
    }
  }
  table <- data.frame(resdf, resdev, c(NA, -diff(resdf)), c(NA, -diff(resdev)))
  variables <- lapply(object, function(x) paste(deparse(x$formula), collapse = "\n"))
  dimnames(table) <- list(1L:nmodels, c("Resid. Df", "Resid. Dev", "Df", "Deviance"))
  if (doscore) table <- cbind(table, Rao = score)
  title <- "Analysis of Deviance Table\n"
  topnote <- paste0("Model ", format(1L:nmodels), ": ", variables, collapse = "\n")
  if(!is.null(test)) {
    bigmodel <- object[[order(resdf)[1L]]]
    dispersion <- summary(bigmodel)$dispersion
    df.dispersion <- if (dispersion == 1)
      Inf
    else min(resdf)
    if (test == "F" && df.dispersion == Inf)
      warning("using F test with a binomial family is inappropriate",
              domain = NA, call. = FALSE)
    table <- stat.anova(table = table, test = test, scale = dispersion,
                        df.scale = df.dispersion, n = length(bigmodel$residuals))
  }
  structure(table, heading = c(title, topnote), class = c("anova", "data.frame"))
}