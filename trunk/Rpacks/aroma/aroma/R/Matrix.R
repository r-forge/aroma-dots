setConstructorS3("Matrix", function(x=0, nrow=1, ncol=1, byrow=FALSE, dimnames=NULL) {
  if (missing(nrow) && missing(ncol) && is.matrix(x)) {
    nrow <- nrow(x);
    ncol <- ncol(x);
  }

  x <- as.vector(x);
  if (length(x) < nrow*ncol)
    x <- rep(x, length.out=nrow*ncol);

  if (missing(nrow))
    nrow <- ceiling(length(x)/ncol)
  else if (missing(ncol))
    ncol <- ceiling(length(x)/nrow);
  
  x <- .Internal(matrix(x, nrow, ncol, byrow));
  x <- as.vector(x);

  extend(MultiwayArray(x, dim=c(nrow,ncol), dimnames=dimnames),
                                                     c("Matrix", "matrix"));
})


setMethodS3("flip", "Matrix", function(this, MARGIN=1) {
  flip <- this;
  if (MARGIN == 1) {
    for (k in 1:ncol(this))
      flip[,k] <- rev(this[,k]);
  } else if (MARGIN == 2) {
    for (k in 1:nrow(this))
      flip[k,] <- rev(this[k,]);
  }
  flip;
})



setMethodS3("rotate", "Matrix", function(this, angle=90) {
  if (diff(dim(this)) != 0)
    stop("Only square matrices can be rotated.");

  if ( (angle %% 90) != 0 )
    stop("Only rotations of 90*n degrees, where n is an integer, is supported.");

  angle <- angle %/% 90;
  angle <- angle %% 4;

  if (angle == 1)
    return(flip(t(this)))
  else if (angle == 2)
    return(flip(t(flip(t(this)))))
  else if (angle == 3)
    return(t(flip(this)));

  this;
})



## setMethodS3("as.DataFrame", "Matrix", function(this, row.names = NULL, optional = FALSE) {
##   DataFrame(as.data.frame(this, row.names=row.names, optional=optional));
## }, deprecated=TRUE)



############################################################################
# HISTORY:
# 2008-01-15
# o CLEAN UP: Removed as.DataFrame() of Matrix.  It was broken anyway.
# 2002-10-24
# o Note that MultiwayArray extends Object with core value. Be careful!!!
# 2002-09-12
# o This class is "identical" to the same in com.braju.sma.
# o Removed the usage of ObjectX2() and extends().
# 2001-11-18
# * Added m o difiers() and generic functions.
# 2001-11-16
# * Added flip() and rotate().
# 2001-11-14
# * Created GSRArray.
# * Created classes MultiwayArray and Matrix.
# * Created!
############################################################################
