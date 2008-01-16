%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{gridSize.Layout}
\alias{gridSize.Layout}
\alias{Layout.gridSize}
\alias{gridSize.Layout}
\alias{gridSize,Layout-method}


\title{Gets the number of spots in each grid}

\description{
 Gets the number of spots in each grid.
}

\usage{\method{gridSize}{Layout}(this, ...)}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
  print(gridSize(layout))   # 324 = 18*18
}
\keyword{internal}
\keyword{methods}
