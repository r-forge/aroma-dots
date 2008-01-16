%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{toXYMatrix.Layout}
\alias{toXYMatrix.Layout}
\alias{Layout.toXYMatrix}
\alias{toXYMatrix.Layout}
\alias{toXYMatrix,Layout-method}


\title{Layouts out values in the same way as the spots are layout}

\description{
 Layouts out values in the same way as the spots are layout.
}

\usage{\method{toXYMatrix}{Layout}(this, values=seq(nbrOfSpots(this)), flip=FALSE, ...)}

\arguments{
  \item{values}{The values to be put into the matrix.}
  \item{flip}{If \code{\link[base:logical]{TRUE}} the resulting matrix is flipped
     vertically before returned.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  layout <- Layout(3,2, 2,3)

  # Print the indices of the each of the spots as they appear
  # on the slide(s).
  print(toXYMatrix(layout))

  #      [,1] [,2] [,3] [,4] [,5] [,6]
  # [1,]    1    2    3    7    8    9
  # [2,]    4    5    6   10   11   12
  # [3,]   13   14   15   19   20   21
  # [4,]   16   17   18   22   23   24
  # [5,]   25   26   27   31   32   33
  # [6,]   28   29   30   34   35   36
}
\keyword{internal}
\keyword{methods}
