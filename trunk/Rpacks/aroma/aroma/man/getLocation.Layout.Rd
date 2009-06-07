%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getLocation.Layout}
\alias{getLocation.Layout}
\alias{Layout.getLocation}
\alias{getLocation.Layout}
\alias{getLocation,Layout-method}


\title{Gets the location of a spot given its index}

\usage{\method{getLocation}{Layout}(this, index=NULL, ...)}

\arguments{
 \item{index}{The spot index of one or many spots to be found. All values
  much be within a valid range otherwise an exception is thrown.}
}

\description{
 Gets the location of a spot given its index. Returns a \code{\link[base]{vector}} containing
 the grid row and the grid column and the spot row and the spot column in
 that grid.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
   layout <- Layout(4,4, 18,18)
   loc <- getLocation(layout, 2001)      # 2 3 4 3
   idx <- getIndex(layout, loc)          # 2001
}

\seealso{
  This method corresponds to \code{\link[sma]{id2image}} in the sma
  package.
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
