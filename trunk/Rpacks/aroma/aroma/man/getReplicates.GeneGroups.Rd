%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GeneGroups.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getReplicates.GeneGroups}
\alias{getReplicates.GeneGroups}
\alias{GeneGroups.getReplicates}
\alias{getReplicates.GeneGroups}
\alias{getReplicates,GeneGroups-method}


\title{Gets all spots that exist with a certain number of replicates}

\usage{\method{getReplicates}{GeneGroups}(this, nbrOfReplicates, value=NULL, asMatrix=TRUE, ...)}

\description{
 Gets all spots that exist with a certain number of replicates.
}

\arguments{
  \item{nbrOfReplicates}{The number of replicates to be selected for, i.e.
    the returned genes all have this number of replicates.}
  \item{value}{A \code{\link[base]{vector}} containing the values to be returned ordered in
    the same way as the spot indicies. If \code{\link[base]{NULL}}, the spot indices are
    returned.}
  \item{asMatrix}{If \code{\link[base:logical]{TRUE}}, the returned values are returned as a
    matrix with \code{nbrOfReplicates} columns and where each row
    represents one gene and the rows are named as the genes,
    otherwise a \code{\link[base]{list}} is returned where each element represents one gene
    and contains \code{nbrOfReplicates} values and the elements are named
    as the genes.}
}

\value{
  Returns a \code{\link[base]{list}} or a \code{\link[base]{matrix}}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  Shorthand versions: \code{\link[aroma:getDuplicates.GeneGroups]{*getDuplicates}()},
  \code{\link[aroma:getTriplicates.GeneGroups]{*getTriplicates}()}, \code{\link[aroma:getQuadruplicates.GeneGroups]{*getQuadruplicates}()}.
  For more information see \code{\link{GeneGroups}}.
}
\keyword{internal}
\keyword{methods}
