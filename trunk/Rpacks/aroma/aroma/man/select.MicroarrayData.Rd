%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{select.MicroarrayData}
\alias{select.MicroarrayData}
\alias{MicroarrayData.select}
\alias{select.MicroarrayData}
\alias{select,MicroarrayData-method}


\title{Selects rows with certain criteria}

\usage{\method{select}{MicroarrayData}(this, ...)}

\arguments{
  \item{...}{The criteria (of fields etc) for returning a row. The indices
             of all matching rows are returned.}
}

\description{
  Given a search criteria, the indices of the rows that fullfills the
  search criteria will be returned.
}

\value{Returns a \code{\link[base]{vector}} of indices.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:extract.MicroarrayData]{*extract}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
