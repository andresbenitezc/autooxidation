####
# Author: Andres Benitez
####
\name{deconSpectra}
\alias{deconSpectra}
\alias{coef.nlsArray}
\alias{resid.nlsArray}
\alias{residuals.nlsArray}
\alias{plot.nlsArray}
\title{Deconvolution of Multiple Spectrums from References}
\description{
The function will deconvolute a given spectrum time series from multiple reference spectrums spectrums provided.
}
\usage{
deconSpectra(S, W, R, M = NULL, baseline = TRUE, lambda = TRUE, ...)
}
\arguments{
  \item{S}{ Matrix with time series of absorbance values for spectrum to be deconvoluted. }
  \item{W}{ Array with values for wavelengths corresponding to absorbance values. }
  \item{R}{ Matrix with the	absorbance values for reference spectrums. }
  \item{M}{ Optional array containing the time values corresponding to the absorbance value time series. }
  \item{baseline}{ Logical indicating if a baseline correction should be used when calculating the spectrum deconvolution. }
  \item{lambda}{ Logical indicating if a lambda function correction should be used when calculating the spectrum deconvolution. }
  \item{\dots}{ Not used. }
}
\value{
  Returns a matrix with relative amounts of each of the reference species present at each time point. If a time series is included that is used returned as the first column of the matrix.
  \item{resid}{Returns residuals for the individual fits.}
  \item{coef }{Returns the coefficients of the fits}
}
\author{ Andres S. Benitez Cardenas }
\examples{
	print("test")
}
\keyword{regression}
