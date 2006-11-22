%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getTags.AffymetrixFile}
\alias{getTags.AffymetrixFile}
\alias{AffymetrixFile.getTags}
\alias{getTags.AffymetrixFile}
\alias{getTags,AffymetrixFile-method}


\title{Gets the tags of the file}

\description{
  Gets the tags of the file.
}

\usage{\method{getTags}{AffymetrixFile}(this, ...)}

\arguments{
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{character}} \code{\link[base]{vector}} or \code{\link[base]{NULL}}.
}

\value{
 The \emph{tags} of a filename are the comma separated parts of the
 filename that follows the the first comma, if any, and that preceeds the
 last period (the filename extension).
 For instance, the tags of \code{path/to/foo,a.2,b.ext} are
 \code{a.2} and \code{b}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:getName.AffymetrixFile]{*getName}()}.
  For more information see \code{\link{AffymetrixFile}}.
}
\keyword{internal}
\keyword{methods}
