%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FirmaModel.R
% 
% on Fri Dec 14 14:17:46 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getFirmaScores.FirmaModel}
\alias{getFirmaScores.FirmaModel}
\alias{FirmaModel.getFirmaScores}
\alias{getFirmaScores.FirmaModel}
\alias{getFirmaScores,FirmaModel-method}


\title{Gets the set of FIRMA results for this model}

\description{
 Gets the set of FIRMA results for this model.
 There is one chip-effect file per array.
}

\usage{\method{getFirmaScores}{FirmaModel}(this, ..., verbose=FALSE)}

\arguments{
  \item{...}{Not used.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link{FirmaSet}} object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{FirmaModel}}.
}
\keyword{internal}
\keyword{methods}
