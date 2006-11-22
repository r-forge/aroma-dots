%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getName.AffymetrixFile}
\alias{getName.AffymetrixFile}
\alias{AffymetrixFile.getName}
\alias{getName.AffymetrixFile}
\alias{getName,AffymetrixFile-method}


\title{Gets the name of the file}

\description{
  Gets the name of the file.
}

\usage{\method{getName}{AffymetrixFile}(this, ...)}

\arguments{
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{character}}.
}

\value{
 The name of a file is the part of the filename without the extension and
 that preceeds any comma.
 For instance, the name of \code{path/to/foo,a.2,b.ext} is \code{foo}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:getTags.AffymetrixFile]{*getTags}()}.
  For more information see \code{\link{AffymetrixFile}}.
}
\keyword{internal}
\keyword{methods}
