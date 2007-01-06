%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.BG.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{bgAdjustRma.AffymetrixCelSet}
\alias{bgAdjustRma.AffymetrixCelSet}
\alias{AffymetrixCelSet.bgAdjustRma}
\alias{bgAdjustRma.AffymetrixCelSet}
\alias{bgAdjustRma,AffymetrixCelSet-method}


\title{Applies RMA background correction to a set of
CEL files}

\description{
 Applies RMA background correction to a set of
CEL files.

 Adapted from \code{\link[affy]{bg.adjust}} in the \pkg{affy} package.
}

\usage{\method{bgAdjustRma}{AffymetrixCelSet}(this, path=NULL, name="bgRma", overwrite=FALSE, skip=!overwrite, ..., verbose=FALSE)}

\arguments{
  \item{path}{The path where to save the adjusted data files.}
  \item{name}{Name of the set containing the background corrected files.}
  \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already adjusted arrays are overwritten,
    unless skipped, otherwise an error is thrown.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the array is not normalized if it already exists.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{...}{Not used.}
}

\value{
 Returns the background adjusted \code{\link{AffymetrixCelFile}} object.
}

\author{
  Ken Simpson (ksimpson[at]wehi.edu.au).
}

\seealso{
 \code{\link[affy]{bg.adjust}}
 For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
