%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ParallelFilter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{changeInput.ParallelFilter}
\alias{changeInput.ParallelFilter}
\alias{ParallelFilter.changeInput}
\alias{changeInput.ParallelFilter}
\alias{changeInput,ParallelFilter-method}


\title{Replaces input filters with other ones}

\usage{\method{changeInput}{ParallelFilter}(this, newInput, oldInput=NULL, recursive=TRUE, ...)}

\arguments{
 \item{recursive}{If \code{\link[base:logical]{TRUE}}, the same method with the arguments are called
   on all input filters.}
}

\description{
 Replaces input filters with other ones.
}

\value{
  Returns nothing.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ParallelFilter}}..
}
\keyword{internal}
\keyword{methods}
