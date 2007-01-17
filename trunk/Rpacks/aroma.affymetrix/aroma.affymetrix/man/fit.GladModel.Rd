%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GladModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fit.GladModel}
\alias{fit.GladModel}
\alias{GladModel.fit}
\alias{fit.GladModel}
\alias{fit,GladModel-method}


\title{Fits the model}

\description{
 Fits the model.
}

\usage{\method{fit}{GladModel}(this, arrays=NULL, chromosomes=getChromosomes(this), force=FALSE, ..., .retResults=FALSE, verbose=FALSE)}

\arguments{
  \item{arrays}{A \code{\link[base]{vector}} of array indices specifying which arrays to
   be considered.  If \code{\link[base]{NULL}}, all are processed.}
  \item{chromosome}{A \code{\link[base]{vector}} of chromosomes indices specifying which
    chromosomes to be considered.  If \code{\link[base]{NULL}}, all are processed.}
  \item{force}{If \code{\link[base:logical]{FALSE}}, the model will not be fitted again if it was
    already fitted.}
  \item{...}{Not used.}
  \item{.retResults}{If \code{\link[base:logical]{TRUE}}, GLAD fit structures are returned for each
    fitted array and chromosome.}
  \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a named \code{\link[base]{list}} of named \code{\link[base]{list}}s.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{GladModel}}.
}
\keyword{internal}
\keyword{methods}
