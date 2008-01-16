%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{nbrOfGrids.Layout}
\alias{nbrOfGrids.Layout}
\alias{Layout.nbrOfGrids}
\alias{nbrOfGrids.Layout}
\alias{nbrOfGrids,Layout-method}


\title{Gets the number of grids on a microarray}

\description{
 Calculates the number of grids on the microarray slide.
}

\usage{\method{nbrOfGrids}{Layout}(this, ...)}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
  print(nbrOfGrids(layout))  # 16
}
\keyword{internal}
\keyword{methods}
