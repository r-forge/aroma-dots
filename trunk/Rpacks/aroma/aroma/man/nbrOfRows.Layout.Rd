%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{nbrOfRows.Layout}
\alias{nbrOfRows.Layout}
\alias{Layout.nbrOfRows}
\alias{nbrOfRows.Layout}
\alias{nbrOfRows,Layout-method}


\title{Gets the number of rows on a microarray}

\description{
 Gets the number of rows on a microarray.
}

\usage{\method{nbrOfRows}{Layout}(this, ...)}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
  print(nbrOfRows(layout))  # 72 = 4*18
}
\keyword{internal}
\keyword{methods}
