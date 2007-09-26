%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SrmaModel.fit.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fit.SrmaModel}
\alias{fit.SrmaModel}
\alias{SrmaModel.fit}
\alias{fit.SrmaModel}
\alias{fit,SrmaModel-method}


\title{Fits the CSRMA model for one chromosome across samples}

\description{
 Fits the CSRMA model for one chromosome across samples.
}

\usage{\method{fit}{SrmaModel}(this, chromosome, force=FALSE, ..., verbose=FALSE)}

\arguments{
  \item{data}{A \code{\link[base]{data.frame}} with columns \code{M} (log-ratio) and
     \code{x} (locus position).
  }
  \item{chromosome}{An \code{\link[base]{integer}} specifying the index of the chromosome to
     be fitted.}
  \item{...}{Additional arguments passed down to the internal fit function.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns the \code{\link[GLAD]{profileCGH}} object returned by \code{\link[GLAD]{glad}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Internally \code{\link[GLAD]{glad}} is used.
  For more information see \code{\link{SrmaModel}}.
}
\keyword{internal}
\keyword{methods}
