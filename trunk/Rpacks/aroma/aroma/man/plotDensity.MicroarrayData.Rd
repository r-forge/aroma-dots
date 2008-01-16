%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.PLOT.R
% 
% on Tue Jan 15 18:33:31 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotDensity.MicroarrayData}
\alias{plotDensity.MicroarrayData}
\alias{MicroarrayData.plotDensity}
\alias{plotDensity.MicroarrayData}
\alias{plotDensity,MicroarrayData-method}


\title{Plots the empirical density distribution of a field}

\description{
 Plots the empirical density distribution of a field.
}

\usage{\method{plotDensity}{MicroarrayData}(this, what, slides=NULL, log=NULL, xlim=NULL, ylim=NULL, xlab=what, ylab="density (integrates to one)", col=NULL, add=FALSE, legend=!add, ...)}

\arguments{
 \item{what}{What to plot. Any field that can be retrieved by
  \code{extract} is accepted.}
 \item{slides}{The slides to be plotted. Each slide will be plotted
  as a separate curve.}
 \item{log}{The log base to be used for taking the logarithm of the
  data before generating the plot. If \code{\link[base]{NULL}} the non-logarithm data
  is plotted.}
 \item{xlim,ylim}{\code{\link[base]{character}} \code{\link[base]{vector}} of length 2. The x and y limits.}
 \item{xlab,ylab}{\code{\link[base]{character}} string for labels on x and y axis.}
 \item{col}{The color(s) of the curves.}
 \item{add}{If \code{\link[base:logical]{TRUE}}, the curves are plotted in the current plot,
  otherwise a new is created.}
 \item{legend}{If \code{\link[base:logical]{TRUE}}, a legend is added, otherwise not.}
 \item{...}{Common arguments accepted by most plot functions.
  For more information see \code{\link[graphics]{par}} and \code{\link[graphics]{plot}}.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:plot.MicroarrayData]{*plot}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
