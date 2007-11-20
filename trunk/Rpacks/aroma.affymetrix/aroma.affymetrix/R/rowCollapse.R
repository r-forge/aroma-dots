###########################################################################/**
# @RdocDefault rowCollapse
#
# @title "Extracts one cell per row from a matrix"
#
# \description{
#  @get "title".
#  The implementation is optimized for memory and speed.
# }
# 
# @synopsis
#
# \arguments{
#   \item{x}{An NxK @matrix.}
#   \item{cols}{An index @vector of (maximum) length N specifying the
#    columns to be extracted.}
#   \item{...}{Not used.}
# }
# 
# \value{
#   @returns a @vector of length N.
# }
#
# \examples{
#   x <- matrix(1:27, ncol=3)
#
#   y <- rowCollapse(x, 1)
#   stopifnot(identical(y, x[,1]))
#
#   y <- rowCollapse(x, 2)
#   stopifnot(identical(y, x[,2]))
#
#   y <- rowCollapse(x, c(1,1,1,1,1,3,3,3,3))
#   stopifnot(identical(y, c(x[1:5,1], x[6:9,3])))
#
#   y <- rowCollapse(x, 1:3)
#   print(y)
#   yT <- c(x[1,1],x[2,2],x[3,3],x[4,1],x[5,2],x[6,3],x[7,1],x[8,2],x[9,3])
#   stopifnot(identical(y, yT))
# }
#
# @author
#
# @keyword utility
#*/########################################################################### 
setMethodS3("rowCollapse", "matrix", function(x, cols, ...) {
  dim <- dim(x);
  colOffsets <- c(0, cumsum(rep(dim[1], dim[2]-1)));
  idxs <- seq_len(dim[1]);
  cols <- rep(cols, length.out=dim[1]);
  colOffsets <- colOffsets[cols];
  idxs <- idxs + colOffsets;
  x[idxs];
})


############################################################################
# HISTORY:
# 2007-10-21
# o Created.
############################################################################
