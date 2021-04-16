logbin.smooth.reparam <- function(coefficients, interpret, type = c("cem","em"), allref, knots,
                                  design.knots, design.param, subset, na.action) {
  type <- match.arg(type)
  coefs.new <- coefficients
  coefs.new.int <- coefs.new[1]
  coefs.new.rm <- NULL
  
  allsmthnames <- NULL
  
  smthterms <- sapply(interpret$smooth.spec, "[[", "term")
  smthnames <- NULL
  for (smth in smthterms) {
    smthlabel <- interpret$smooth.spec[[smth]]$termlabel
    which.smth <- which(substr(names(coefficients), 1, nchar(smthlabel)) == smthlabel)
    coefs.smth <- -coefficients[which.smth]
    smthnames <- c(smthnames, names(coefs.smth))
    smthtype <- class(interpret$smooth.spec[[smth]])
    if (smthtype == "Iso.smooth") {
      coefs.new.int <- coefs.new.int + sum(coefs.smth)
      coefs.smth.new <- -coefs.smth
      names.smth.new <- names(coefs.smth)
    } else if (smthtype == "B.smooth") {
      ref <- allref$allref[[smth]][[as.numeric(design.param[smth])]]
      num.knots <- length(knots[[smth]])
      if (length(ref) == 1) {
        if (type == "cem") {
          coefs.smth.temp <- append(coefs.smth, 0, after = ref - 1)
          coefs.new.int <- coefs.new.int + coefs.smth.temp[1]
          coefs.smth.new <- (coefs.smth.temp - coefs.smth.temp[1])[-1]
        } else {
          coefs.smth.temp <- coefs.smth
          coefs.new.int <- coefs.new.int + coefs.smth.temp[1]
          coefs.smth.new <- (coefs.smth.temp - coefs.smth.temp[1])[-1]
        }
      } else {
        coefs.smth.temp <- c(rev(cumsum(rev(coefs.smth))), 0)
        coefs.smth.ord <- coefs.smth.temp[order(ref)]
        coefs.new.int <- coefs.new.int + coefs.smth.ord[1]
        coefs.smth.new <- (coefs.smth.ord - coefs.smth.ord[1])[-1]
      }
      names.smth.new <- paste0(smthlabel, 2:(num.knots-3))
    }
    if (smthtype == "B.smooth" && length(ref) == 1 && type == "em") {
      coefs.new[which.smth] <- c(0,coefs.smth.new)
      names(coefs.new)[which.smth] <- c("TEMP",names.smth.new)
      coefs.new.rm <- c(coefs.new.rm, which.smth[1])
    } else {
      coefs.new[which.smth] <- coefs.smth.new
      names(coefs.new)[which.smth] <- names.smth.new
    }
    allsmthnames <- c(allsmthnames, names.smth.new)
  }
  
  coefs.new[1] <- coefs.new.int
  if (!is.null(coefs.new.rm)) coefs.new <- coefs.new[-coefs.new.rm]
  
  dummy.design <- rep(1, length(design.param))
  names(dummy.design) <- names(design.param)
  dummy.allref <- allref
  for (smth in names(allref$allref))
    dummy.allref$allref[[smth]][[1]] <- 1
  modelspec <- logbin.smooth.design(interpret, "cem", dummy.allref, design.knots, dummy.design)
  data.new <- modelspec$data
  data.new[, allsmthnames] <- -data.new[, allsmthnames]
  dummy.frame.call <- call("model.frame", formula = eval(modelspec$formula), data = as.name("data.new"))
  dummy.frame.call$drop.unused.levels <- TRUE
  if (!missing(na.action)) dummy.frame.call$na.action <- na.action
  if (!missing(subset)) dummy.frame.call$subset = subset
  dummy.frame <- eval(dummy.frame.call)
  dummy.terms <- attr(dummy.frame, "terms")
  
  design <- model.matrix(dummy.terms, dummy.frame)
  
  list(coefs = coefs.new, mf = dummy.frame, design = design, mt = dummy.terms,
       smoothnames = smthnames)
}