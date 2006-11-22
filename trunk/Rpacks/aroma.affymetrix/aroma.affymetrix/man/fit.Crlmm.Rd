%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Crlmm.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fit.Crlmm}
\alias{fit.Crlmm}
\alias{Crlmm.fit}
\alias{fit.Crlmm}
\alias{fit,Crlmm-method}


\title{Fits the model}

\description{
 Fits the model.
}

\usage{\method{fit}{Crlmm}(this, recalibrate=TRUE, ..., verbose=FALSE)}

\arguments{
  \item{recalibrate}{If \code{\link[base:logical]{TRUE}}, a second round of the fitting is done.}
  \item{...}{Additional arguments.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[oligo]{SnpCallSetPlus-class}} object.
}

\details{
  Simple benchmarking on shadowfax:
  For 90 CEPH Xba chips it takes ~90 minutes, i.e. 60 s/array.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[oligo]{crlmm}}.
  For more information see \code{\link{Crlmm}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
