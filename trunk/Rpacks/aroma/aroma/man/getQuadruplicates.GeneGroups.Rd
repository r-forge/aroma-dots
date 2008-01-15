%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GeneGroups.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getQuadruplicates.GeneGroups}
\alias{getQuadruplicates.GeneGroups}
\alias{GeneGroups.getQuadruplicates}
\alias{getQuadruplicates.GeneGroups}
\alias{getQuadruplicates,GeneGroups-method}


\title{Gets all quadruplicated spots}

\usage{\method{getQuadruplicates}{GeneGroups}(this, ...)}

\description{
 Gets all quadruplicated spots. This is a shortcut for calling
 \code{\link[aroma:getReplicates.GeneGroups]{*getReplicates}()} with \code{nbrOfReplicates=4}.
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
