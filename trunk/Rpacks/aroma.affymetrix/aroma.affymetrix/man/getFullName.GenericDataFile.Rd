%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenericDataFile.R
% 
% on Mon Mar 31 17:07:46 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getFullName.GenericDataFile}
\alias{getFullName.GenericDataFile}
\alias{GenericDataFile.getFullName}
\alias{getFullName.GenericDataFile}
\alias{getFullName,GenericDataFile-method}


\title{Gets the full name of the file}

\description{
  Gets the full name of the file, that is the filename without the extension.
}

\usage{\method{getFullName}{GenericDataFile}(this, aliased=FALSE, ...)}

\arguments{
 \item{aliased}{If \code{\link[base:logical]{TRUE}}, and an alias has been set, the alias is
    returned, otherwise the default full name is returned.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{character}}.
}

\details{
 The full name of a file is the filename excluding any
 extension (and period).
 For instance, the full name of \code{path/to/foo,a.2,b.ext} is
 \code{foo,a.2,b}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:getName.GenericDataFile]{*getName}()}.
  \code{\link[aroma.affymetrix:getTags.GenericDataFile]{*getTags}()}.
  For more information see \code{\link{GenericDataFile}}.
}
\keyword{internal}
\keyword{methods}
