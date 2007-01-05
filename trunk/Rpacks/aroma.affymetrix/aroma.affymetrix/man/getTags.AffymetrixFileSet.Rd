%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% on Sat Jan 06 01:47:51 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getTags.AffymetrixFileSet}
\alias{getTags.AffymetrixFileSet}
\alias{AffymetrixFileSet.getTags}
\alias{getTags.AffymetrixFileSet}
\alias{getTags,AffymetrixFileSet-method}


\title{Gets the tags of the file set}

\description{
  Gets the tags of the file set.
  Any tag equals \code{"*"} is replaced by the comma separated tags part of
  the file-set pathname.
}

\usage{\method{getTags}{AffymetrixFileSet}(this, ...)}

\arguments{
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{character}} \code{\link[base]{vector}} or \code{\link[base]{NULL}}.
}

\details{
 The \emph{tags} of a file set are the comma separated parts of the
 filename that follows the the first comma, if any, and that preceeds the
 last period (the filename extension).
 For instance, the tags of \code{path/to/foo,a.2,b.ext} are
 \code{a.2} and \code{b}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:setTags.AffymetrixFileSet]{*setTags}()}.
  \code{\link[aroma.affymetrix:getName.AffymetrixFileSet]{*getName}()}.
  For more information see \code{\link{AffymetrixFileSet}}.
}
\keyword{internal}
\keyword{methods}
