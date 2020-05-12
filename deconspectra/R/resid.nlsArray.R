####
# Creates functions for fitting residuals if needed for deconspectra.
# Author: Andres Benitez
####

resid.nlsArray <- function(x) { sapply(x$fits, residuals) }
residuals.nlsArray <- function(x) { sapply(x$fits, residuals) }
