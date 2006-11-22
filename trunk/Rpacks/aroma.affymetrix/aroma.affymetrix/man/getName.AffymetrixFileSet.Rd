%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixFileSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getName.AffymetrixFileSet}
\alias{getName.AffymetrixFileSet}
\alias{AffymetrixFileSet.getName}
\alias{getName.AffymetrixFileSet}
\alias{getName,AffymetrixFileSet-method}


\title{Gets the name of the file set}

\description{
  Gets the name of the file set.
}

\usage{\method{getName}{AffymetrixFileSet}(this, ...)}

\arguments{
 \item{...}{Additional arguments passed to \code{\link[aroma.affymetrix:getFullName.AffymetrixFileSet]{*getFullName}()}, if
   the name is inferred from the pathname.}
}

\value{
  Returns a \code{\link[base]{character}}.
}

\details{
 If the name has not been set explicitly, the name is inferred from the
 pathname of the file set.
 The \emph{name} of a file set is the part of the directory name that
 preceeds the first comma, if any.
 For instance, the name of the file set named \code{foo,a,b} is \code{foo}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:setName.AffymetrixFileSet]{*setName}()}.
  \code{\link[aroma.affymetrix:getFullName.AffymetrixFileSet]{*getFullName}()}.
  \code{\link[aroma.affymetrix:getBranches.AffymetrixFileSet]{*getBranches}()}.
  For more information see \code{\link{AffymetrixFileSet}}.
}
\keyword{internal}
\keyword{methods}
