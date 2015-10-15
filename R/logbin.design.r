logbin.design <- function(terms, data, type = c("cem","em"), allref, mono, design.ref) {
  type <- match.arg(type)
  terms.new <- terms
	data.new <- data
	termlabels <- attr(terms, "term.labels")
  design.type <- sapply(allref, attr, "type")
  ref.vector <- as.vector(design.ref, mode = "integer")
  for (i in seq_len(length(design.ref))) {
    varname <- gsub("`", "", termlabels[i])
    varref <- allref[[termlabels[i]]][[ref.vector[i]]]
    if (design.type[i] == 1) {
      cont.min <- min(data[[varname]], na.rm = TRUE)
      cont.max <- max(data[[varname]], na.rm = TRUE)
      if (varref == 1) data.new[[varname]] <- data[[varname]] - cont.min
      else data.new[[varname]] <- cont.max - data[[varname]]
      if (type == "em" && !mono[termlabels[i]]) {
        varname.rev <- paste0(varname, ".rev")
        if (varref == 1) data.new[[varname.rev]] <- cont.max - data[[varname]]
        else data.new[[varname.rev]] <- data[[varname]] - cont.min
        terms.new <- terms(update(terms.new, as.formula(paste(".~. +", varname.rev))))
      }
    } else if(design.type[i] == 2) {
      data.new[[varname]] <- relevel(factor(data[[varname]]), ref = varref)
      if(type == "cem")
        contrasts(data.new[[varname]]) <- contr.treatment(levels(data.new[[varname]]), base = 1)
      else
        contrasts(data.new[[varname]], nlevels(data.new[[varname]])) <- contr.treatment(levels(data.new[[varname]]), base = 1, contrasts = FALSE)
    } else if(design.type[i] == 3) {
      data.new[[varname]] <- factor(data[[varname]])
      contrasts(data.new[[varname]]) <- contr.isotonic.rev(levels(data.new[[varname]]), perm = varref)
    }
  }
  attr(data.new, "terms") <- terms.new
  X <- model.matrix(terms.new, data.new)	
	X
}