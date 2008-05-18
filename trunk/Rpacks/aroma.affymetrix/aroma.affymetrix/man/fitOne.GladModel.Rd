%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GladModel.fitOne.R
% 
% on Sat May 17 11:54:49 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fitOne.GladModel}
\alias{fitOne.GladModel}
\alias{GladModel.fitOne}
\alias{fitOne.GladModel}
\alias{fitOne,GladModel-method}


\title{Fits the GLAD model for one chromosome in one sample}

\description{
 Fits the GLAD model for one chromosome in one sample.
}

\usage{\method{fitOne}{GladModel}(this, data, chromosome, ..., verbose=FALSE)}

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
  For more information see \code{\link{GladModel}}.
}
\keyword{internal}
\keyword{methods}
