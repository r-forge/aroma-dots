%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  gcrmaSummary.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{gcrmaSummary.AffymetrixCelSet}
\alias{gcrmaSummary.AffymetrixCelSet}
\alias{AffymetrixCelSet.gcrmaSummary}
\alias{gcrmaSummary.AffymetrixCelSet}
\alias{gcrmaSummary,AffymetrixCelSet-method}


\title{Calculates the GCRMA expression summary}

\description{
 Calculates the GCRMA expression summary.

Applies GC-based background correction, followed by quantile normalisation
and then a summary by either median polish or robust regression.
}

\usage{\method{gcrmaSummary}{AffymetrixCelSet}(this, path=NULL, name="gcrma", bgPath=NULL, normPath=NULL, summaryMethod="rlm", probePath=NULL, affinities=NULL, type="fullmodel", indicesNegativeControl=NULL, opticalAdjust=TRUE, gsbAdjust=TRUE, k=6 * fast + 0.5 * (1 - fast), rho=0.7, stretch=1.15 * fast + (1 - fast), fast=TRUE, ..., verbose=FALSE)}

\arguments{
  \item{path}{The path to save the expression summaries.}
  \item{name}{Name of the data set containing the expression summaries.}
  \item{bgPath}{Directory in which to store background-adjusted signals.}
  \item{normPath}{Directory in which to store normalised signals.}
  \item{summaryMethod}{Either "medianpolish" or "rlm".}
  \item{type}{The type of background correction.  Currently accepted types
      are "fullmodel" (the default, uses MMs) and "affinities" (uses
      probe sequence only).}
  \item{indicesNegativeControl}{Locations of any negative control
      probes (e.g., the anti-genomic controls on the human exon array).
      If \code{\link[base]{NULL}} and type=="affinities", MMs are used as the negative
      controls.}
  \item{opticalAdjust}{If \code{\link[base:logical]{TRUE}}, apply correction for optical effect,
      as in \code{\link[gcrma]{bg.adjust.optical}}.}
  \item{gsbAdjust}{Should we adjust for specific binding (defaults to
       \code{\link[base:logical]{TRUE}})?}
  \item{k}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
  \item{rho}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
  \item{stretch}{Tuning parameter passed to \code{gcrma::bg.adjust.gcrma}.}
  \item{fast}{If \code{\link[base:logical]{TRUE}}, an ad hoc transformation of the PM is performed
      (\code{gcrma::gcrma.bg.transformation.fast}).}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already adjusted arrays are overwritten,
    unless skipped, otherwise an error is thrown.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the array is not normalized if it already exists.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link{AffymetrixCelSet}} containing the expression summaries.
}

\author{
  Ken Simpson (ksimpson[at]wehi.edu.au).
}

\seealso{
 \code{\link[gcrma]{bg.adjust.gcrma}}
 For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
