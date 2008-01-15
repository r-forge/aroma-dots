%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffineModelFit.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AffineModelFit}
\docType{class}
\alias{AffineModelFit}


\title{Model fit of affine microarray models}

\description{
  Model fit of affine microarray models.
}

\usage{AffineModelFit(a, b, adiag, fitted=NULL, weights=NULL, residuals=NULL, eigen=NULL, y=NULL, converged=NULL, nbrOfIterations=NULL, t0=NULL, t=NULL, ...)}

\arguments{
 \item{a,b,adiag}{Parameters.}
 \item{residuals}{.}
 \item{eigen}{.}
 \item{converged}{If \code{\link[base:logical]{TRUE}}, it indicates that the estimate converged. Otherwise, it did not.}
 \item{nbrOfIterations}{Number of iterations until convergence or stop.}
 \item{t,t0}{.}
 \item{fitted}{a \code{\link[base]{vector}} of fitted values.}
 \item{weights}{(only for weighted fits) a \code{\link[base]{vector}} the weights used for
   the fit. If a iterative re-weighted method is used, these are the
   weights used in the last iteration.}
 \item{y}{if requested, the response \code{\link[base]{matrix}} used.}
 \item{...}{Other named arguments stored.}
}

\value{
  Returns a \code{\link[base]{list}} containing the fitted model.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
}

\seealso{
  \code{\link{RGData.calibrateMultiscan}}
  \code{\link{RGData.normalizeAffine}}
  For more information see \code{\link{AffineModelFit}}.
}
\keyword{classes}
