%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getPosition.Layout}
\alias{getPosition.Layout}
\alias{Layout.getPosition}
\alias{getPosition.Layout}
\alias{getPosition,Layout-method}


\title{Gets the position of a set of spots given their indices}

\usage{\method{getPosition}{Layout}(this, index=NULL, ...)}

\arguments{
 \item{index}{The spot index of one or many spots to be found. All values
  much be within a valid range otherwise an exception is thrown.
  If \code{\link[base]{NULL}}, all spots are considered.}
}

\description{
 Gets the position of a set of spots given their indices. The position of
 a spot is the pair (row, column) where row is the row in the microarray
 and column is the column in the microarray where the spot is positioned.
}

\value{
 Returns a "named" \code{\link[base]{matrix}} where each row contains the row and the
 column of the spots.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
   layout <- Layout(4,4, 18,18)
   xy <- getPosition(layout, 2001)           # 22 39
   xy <- getPosition(layout, c(2001,2002))   # 22 39; 22 40
}

\seealso{
  \code{\link[aroma:getIndex.Layout]{*getIndex}()}, \code{\link[aroma:getIndices.Layout]{*getIndices}()},
  \code{\link[aroma:getLocation.Layout]{*getLocation}()}.
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
