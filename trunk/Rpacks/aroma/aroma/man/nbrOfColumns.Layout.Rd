%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{nbrOfColumns.Layout}
\alias{nbrOfColumns.Layout}
\alias{Layout.nbrOfColumns}
\alias{nbrOfColumns.Layout}
\alias{nbrOfColumns,Layout-method}


\title{Gets the number of columns on a microarray}

\description{
 Gets the number of columns on a microarray.
}

\usage{\method{nbrOfColumns}{Layout}(this, ...)}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  layout <- Layout(ngrid.r=4, ngrid.c=4, nspot.r=18, nspot.c=18)
  print(nbrOfColumns(layout))  # 72 = 4*18
}
\keyword{internal}
\keyword{methods}
