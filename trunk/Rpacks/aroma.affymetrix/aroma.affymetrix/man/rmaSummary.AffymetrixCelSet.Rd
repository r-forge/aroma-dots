%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  rmaSummary.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{rmaSummary.AffymetrixCelSet}
\alias{rmaSummary.AffymetrixCelSet}
\alias{AffymetrixCelSet.rmaSummary}
\alias{rmaSummary.AffymetrixCelSet}
\alias{rmaSummary,AffymetrixCelSet-method}


\title{Calculates the RMA expression summary}

\description{
 Calculates the RMA expression summary.

Models the observed signal as normal background + exponential signal, then
applies quantile normalisation followed by a summary by either median
polish or robust regression.

}

\usage{\method{rmaSummary}{AffymetrixCelSet}(this, path=NULL, name="rma", bgPath=NULL, normPath=NULL, summaryMethod="rlm", ..., verbose=FALSE)}

\arguments{
  \item{path}{The path to save the expression summaries.}
  \item{name}{Name of the data set containing the expression summaries.}
  \item{bgPath}{Directory in which to store background-adjusted signals.}
  \item{normPath}{Directory in which to store normalised signals.}
  \item{summaryMethod}{Either "medianpolish" or "rlm".}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already adjusted arrays are overwritten,
    unless skipped, otherwise an error is thrown.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the array is not normalized if it already exists.}
  \item{...}{Additional arguments passed to some of the internal functions.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link{AffymetrixCelSet}} containing the expression summaries.
}

\author{
  Ken Simpson (ksimpson[at]wehi.edu.au).
}

\seealso{
 See package \pkg{affy}.
 For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
