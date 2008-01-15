############################################################################
# Private utility functions for com.braju.sma
#
# Author(s): Henrik Bengtsson, henrikb@braju.com
# Date     : March 2001
############################################################################
# Gets the dimension of a matrix as a string of format "(rows,cols)".
com.braju.sma.dimStr <- function(matrix) {
  paste(sep="", "(", nrow(matrix), "x", ncol(matrix), ")")
}


setMethodS3("matrixToList", "default", function(mat, na.rm=FALSE) {
  l0 <- apply(mat, MARGIN=1, FUN=function(row) list(if (na.rm) na.omit(row) else row))
  l1 <- lapply(l0, FUN=function(row) unlist(row))
})


setMethodS3("listToMatrix", "default", function(list, fill.value=NA) {
  ncol <- max(unlist(lapply(list, FUN=length)));
  filled.list <- lapply(list, FUN = function(row) {
                                      len <- length(row);
                                      c(row, rep(fill.value, ncol-len));
                                    });
  matrix(unlist(filled.list), ncol=ncol, byrow=TRUE);
})



setMethodS3("is.above", "default", function(line, xy, na.replace=FALSE) {
  x.idx <- apply(as.matrix(xy$x), MARGIN=1, FUN=function(x) {
      x0 <- which(x >= line$x);
      if (length(x0) == 0) 1 else max(x0);
    }
  )

  xs <- matrix(c(line$x[x.idx], line$x[x.idx+1]), nrow=2, byrow=TRUE); 
  delta <- (xy$x-xs[1,])/(xs[2,]-xs[1,]);

  ys <- matrix(c(line$y[x.idx], line$y[x.idx+1]), nrow=2, byrow=TRUE);

  y.at.x <- (ys[2,]-ys[1,])*delta+ys[1,];
  y.dist <- xy$y-y.at.x;

  res <- (y.dist > 0);
  if (length(res) > 0 && !is.null(na.replace)) {
    res[is.na(res)] <- na.replace;
  }
  res;
})


############################################################################
# HISTORY:
# 2003-04-25
# o Updated with setMethodS3().
# 2001-07-15
# * Added matrixToList and listToMatrix. Sooner to be moved to R.basic.
# 2001-03-25
# * Created.
############################################################################


