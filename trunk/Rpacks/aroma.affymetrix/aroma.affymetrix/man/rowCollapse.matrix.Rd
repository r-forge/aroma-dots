%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  rowCollapse.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{rowCollapse.matrix}
\alias{rowCollapse.matrix}
\alias{rowCollapse.matrix}


 \title{Extracts one cell per row from a matrix}

 \description{
  Extracts one cell per row from a matrix.
  The implementation is optimized for memory and speed.
 }

 \usage{\method{rowCollapse}{matrix}(x, cols, ...)}

 \arguments{
   \item{x}{An NxK \code{\link[base]{matrix}}.}
   \item{cols}{An index \code{\link[base]{vector}} of (maximum) length N specifying the
    columns to be extracted.}
   \item{...}{Not used.}
 }

 \value{
   @returns a \code{\link[base]{vector}} of length N.
 }

 \examples{
   x <- matrix(1:27, ncol=3)

   y <- rowCollapse(x, 1)
   stopifnot(identical(y, x[,1]))

   y <- rowCollapse(x, 2)
   stopifnot(identical(y, x[,2]))

   y <- rowCollapse(x, c(1,1,1,1,1,3,3,3,3))
   stopifnot(identical(y, c(x[1:5,1], x[6:9,3])))

   y <- rowCollapse(x, 1:3)
   print(y)
   yT <- c(x[1,1],x[2,2],x[3,3],x[4,1],x[5,2],x[6,3],x[7,1],x[8,2],x[9,3])
   stopifnot(identical(y, yT))
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{methods}
\keyword{utilities}
