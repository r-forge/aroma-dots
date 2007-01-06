%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CrlmmModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{classifyAsXorXX.CrlmmModel}
\alias{classifyAsXorXX.CrlmmModel}
\alias{CrlmmModel.classifyAsXorXX}
\alias{classifyAsXorXX.CrlmmModel}
\alias{classifyAsXorXX,CrlmmModel-method}


\title{Classifies the samples to have one or two copies of chromosome X}

\description{
 Classifies the samples to have one or two copies of chromosome X.
}

\usage{\method{classifyAsXorXX}{CrlmmModel}(this, ...)}

\arguments{
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{vector}} of length N, where N is the number of samples.
}

\details{
 This method is adopted from the internal \code{snpGenderCall()} function
 in the \pkg{oligo} package.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{CrlmmModel}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
