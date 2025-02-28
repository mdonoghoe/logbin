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

anova.logbin <- function(object, ..., test = NULL)
{
    dotargs <- list(...)
    named <- if (is.null(names(dotargs))) rep(FALSE, length(dotargs))
                else (names(dotargs) != "")
    if (any(named))
        warning("the following arguments to 'anova.logbin' are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
    dotargs <- dotargs[!named]
    is.logbin <- unlist(lapply(dotargs, function(x) inherits(x, "logbin")))
    dotargs <- dotargs[is.logbin]
    if (length(dotargs))
        return(anova.logbinlist(c(list(object), dotargs), test = test))
    else
        stop('anova.logbin does not support anova for a single model. Fit nested models manually and input to anova.logbin')
}