%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Model.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getPath.Model}
\alias{getPath.Model}
\alias{Model.getPath}
\alias{getPath.Model}
\alias{getPath,Model-method}


\title{Gets the path of this model}

\description{
 Gets the path of this model where the parameter files are located.
}

\usage{\method{getPath}{Model}(this, mkdirs=TRUE, ...)}

\arguments{
 \item{mkdirs}{If \code{\link[base:logical]{TRUE}}, the directory is created, if missing.}
 \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{character}} string.
}

\details{
 If the path does not exist, it is created.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{Model}}.
}
\keyword{internal}
\keyword{methods}
