%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{size.Layout}
\alias{size.Layout}
\alias{Layout.size}
\alias{size.Layout}
\alias{size,Layout-method}


\title{Gets the size of a microarray}

\description{
 Calculates the total number of spots on the microarray slide.
}

\usage{\method{size}{Layout}(this, ...)}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
  print(size(layout))  # 5184
}

\seealso{
  \code{\link[aroma:nbrOfSpots.Layout]{*nbrOfSpots}()}.
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
