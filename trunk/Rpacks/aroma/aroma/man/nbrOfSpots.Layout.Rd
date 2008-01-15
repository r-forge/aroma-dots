%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{nbrOfSpots.Layout}
\alias{nbrOfSpots.Layout}
\alias{Layout.nbrOfSpots}
\alias{nbrOfSpots.Layout}
\alias{nbrOfSpots,Layout-method}


\title{Gets the size of a microarray}

\description{
 Calculates the total number of spots on the microarray slide.
}

\usage{\method{nbrOfSpots}{Layout}(this, ...)}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
  print(nbrOfSpots(layout))  # 5184
}
\keyword{internal}
\keyword{methods}
