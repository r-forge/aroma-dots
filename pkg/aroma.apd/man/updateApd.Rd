%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  updateApd.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{updateApd}
\alias{updateApd.default}
\alias{updateApd}

\title{Updates an Affymetrix probe data (APD) file}

\usage{
\method{updateApd}{default}(filename, indices=NULL, data, writeMap=NULL, ..., verbose=FALSE,
  .checkArgs=TRUE)
}

\description{
  Updates an Affymetrix probe data (APD) file.
}

\arguments{
  \item{filename}{The filename of the APD file.}
  \item{indices}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of cell (probe) indices specifying
    which cells to updated.}
  \item{data}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of data elements to be assigned.}
  \item{writeMap}{A \code{\link[base]{vector}} of indicies used to change the order how
     data elements are \emph{written} (by default).}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{.checkArgs}{If \code{\link[base:logical]{TRUE}}, arguments are checked, otherwise not.}
}

\value{
  Returns (invisibly) the pathname of the file updated.
}

\author{Henrik Bengtsson}

\examples{\dontrun{#See ?createApd for an example.}}

\seealso{
  \code{\link{createApd}}() and \code{\link{updateApd}}().
}



\keyword{file}
\keyword{IO}
