%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% on Sat Jan 06 01:47:51 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getFullName.AffymetrixFileSet}
\alias{getFullName.AffymetrixFileSet}
\alias{AffymetrixFileSet.getFullName}
\alias{getFullName.AffymetrixFileSet}
\alias{getFullName,AffymetrixFileSet-method}


\title{Gets the full name of the file set}

\description{
  Gets the full name of the file set, that is the name of the directory without parent directories.
}

\usage{\method{getFullName}{AffymetrixFileSet}(this, parent=1, ...)}

\arguments{
 \item{parent}{The number of generations up in the directory tree the
   directory name should be retrieved.  By default the current directory
   is used.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{character}}.
}

\details{
 By default, the full name of a file set is the name of the directory
 containing all the files, e.g. the name of file set \code{path/to,a,b/*}
 is \code{to,a,b}.
 Argument \code{parent=1} specifies that the parent directory should be
 used, and so on.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixFileSet}}.
}
\keyword{internal}
\keyword{methods}
