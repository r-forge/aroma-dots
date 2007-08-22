%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GladModel.fitOne.R
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
  \item{data}{}
  \item{...}{}
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
