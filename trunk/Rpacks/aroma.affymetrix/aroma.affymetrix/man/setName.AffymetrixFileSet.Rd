%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setName.AffymetrixFileSet}
\alias{setName.AffymetrixFileSet}
\alias{AffymetrixFileSet.setName}
\alias{setName.AffymetrixFileSet}
\alias{setName,AffymetrixFileSet-method}


\title{Sets the name of the file set}

\description{
  Sets the name of the file set.
}

\usage{\method{setName}{AffymetrixFileSet}(this, name=NULL, ...)}

\arguments{
 \item{name}{A \code{\link[base]{character}} string for the new name of the file set.
  The name must consists of valid filename characters, and must not
  contain commas, which are used to separate tags.
  If \code{\link[base]{NULL}}, the name will be inferred from the underlying pathname.}
 \item{...}{Not used.}
}

\value{
  Returns nothing.
}

\details{
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:getName.AffymetrixFileSet]{*getName}()}.
  For more information see \code{\link{AffymetrixFileSet}}.
}
\keyword{internal}
\keyword{methods}
