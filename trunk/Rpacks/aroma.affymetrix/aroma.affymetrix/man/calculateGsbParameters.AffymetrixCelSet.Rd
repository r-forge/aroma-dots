%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.BG.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{calculateGsbParameters.AffymetrixCelSet}
\alias{calculateGsbParameters.AffymetrixCelSet}
\alias{AffymetrixCelSet.calculateGsbParameters}
\alias{calculateGsbParameters.AffymetrixCelSet}
\alias{calculateGsbParameters,AffymetrixCelSet-method}


\title{Computes parameters for adjustment of specific binding}

\description{
 Computes parameters for adjustment of specific binding.
}

\usage{\method{calculateGsbParameters}{AffymetrixCelSet}(this, nbrOfPms=25000, affinities=NULL, path=NULL, ..., verbose=FALSE)}

\arguments{
  \item{nbrOfPms}{The number of random PMs to use in estimation.}
  \item{affinities}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of probe affinities.}
  \item{path}{If an affinities vector is not specified,
     gives the path to a file storing the affinities.}
}

\value{
}

\author{
  Ken Simpson (ksimpson[at]wehi.edu.au).
}
\keyword{internal}
\keyword{methods}
