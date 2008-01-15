%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GeneGroups.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getDuplicates.GeneGroups}
\alias{getDuplicates.GeneGroups}
\alias{GeneGroups.getDuplicates}
\alias{getDuplicates.GeneGroups}
\alias{getDuplicates,GeneGroups-method}


\title{Gets all duplicated spots}

\usage{\method{getDuplicates}{GeneGroups}(this, ...)}

\description{
 Gets all duplicated spots. This is a shortcut for calling
 \code{\link[aroma:getReplicates.GeneGroups]{*getReplicates}()} with \code{nbrOfReplicates=2}.
}

\arguments{
  \item{...}{Arguments accepted by \code{\link[aroma:getReplicates.GeneGroups]{*getReplicates}()}.}
}

\value{
  Returns a \code{\link[base]{list}} or a \code{\link[base]{matrix}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:getReplicates.GeneGroups]{*getReplicates}()}.
  For more information see \code{\link{GeneGroups}}.
}
\keyword{internal}
\keyword{methods}
