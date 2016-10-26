plot.logbin.smooth <- function(x, type = c("response", "link", "diagnostics"), at = data.frame(), 
                        knotlines = TRUE, nobs = 1000, ...) {
  type <- match.arg(type)
  extraargs <- list(...)
  extraargs2 <- extraargs
  extraargs2$ylim <- NULL
  extraargs2$ylab <- NULL
  extraargs2$col <- NULL
  if (type == "diagnostics") plot(structure(x, class = c("glm", "lm")), ...)
  else {
    gp <- interpret.logbin.smooth(x$full.formula)
    allvars <- names(get_all_vars(delete.response(gp$terms), data = x$data))
    smthterms <- names(gp$smooth.spec)
    for (i in 1:length(smthterms)) {
      allvars2 <- allvars[allvars != smthterms[i]]
      smthtype <- class(gp$smooth.spec[[smthterms[i]]])
      ltype = "l"
      if (smthtype == "Iso.smooth")
        ltype = "S"
      smthknots <- x$knots[[smthterms[i]]]
      smthx <- seq(smthknots[1], smthknots[length(smthknots)], len = nobs)
      dat <- list()
      at.new <- at
      at.new[[smthterms[i]]] <- NULL
      if (!setequal(intersect(allvars2, names(at.new)), allvars2))
        stop(gettextf("error predicting for %s: 'at' must contain %s",
                      smthterms[i], paste(allvars2, collapse = ", ")), domain = NA)
      if (nrow(at.new) == 0)           
        dat[[1]] <- setNames(data.frame(tmp = smthx), smthterms[i])
      else {
        for (j in 1:nrow(at.new))
          dat[[j]] <- cbind(setNames(data.frame(tmp = smthx), smthterms[i]),
                            setNames(data.frame(tmp = at.new[j,]), names(at.new)), row.names = NULL)
      }
      pred <- list()
      for (j in 1:length(dat))
        pred[[j]] <- predict(x, newdata = dat[[j]], type = type)
      ylim2 <- c(min(sapply(pred, min)), max(sapply(pred, max)))
      ylab2 <- type
      plotargs <- list(x = smthx, y = pred[[1]], type = "n", xlab = smthterms[i])
      plotargs$ylim <- if(!is.null(extraargs[["ylim"]])) extraargs[["ylim"]] else ylim2
      plotargs$ylab <- if(!is.null(extraargs[["ylab"]])) extraargs[["ylab"]] else ylab2
      do.call(plot, c(plotargs, extraargs2))
      if(knotlines & smthtype != "Iso.smooth")
        abline(v = smthknots, col = "gray", lwd = 0.5)
      plotcols <- rainbow(length(pred))
      if(!is.null(extraargs[["col"]])) plotcols <- rep(extraargs[["col"]], length(pred))
      lineargs <- list(x = smthx, type = ltype)
      for (j in 1:length(pred)) {
        lineargs$y <- pred[[j]]
        lineargs$col <- plotcols[j]
        do.call(lines, c(lineargs, extraargs2))
      }
    }
  }
}