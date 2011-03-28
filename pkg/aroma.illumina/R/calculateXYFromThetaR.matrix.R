setMethodS3("calculateThetaRFromXY", "array", function(xy, ...) {
  # Argument 'xy':
  dim <- dim(xy);
  stopifnot(length(dim) == 3);
  stopifnot(dim[2] == 2);

  # (X,Y) -> (theta,R)
  X <- xy[,1,];
  Y <- xy[,2,];
  R <- X+Y;
  theta <- (2/pi)*atan(Y/X);

  # Allocate results
  res <- xy;
  res[,1,] <- theta;
  res[,2,] <- R;
  dimnames(res)[[2]] <- c("theta", "R");

  res;
}) # calculateThetaRFromXY()


setMethodS3("calculateThetaRFromXY", "matrix", function(xy, ...) {
  # Argument 'xy':
  dim <- dim(xy);
  dim(xy) <- c(dim, 1L);
  res <- calculateThetaRFromXY(xy, ...);
  dim(res) <- dim;
  res;
})


############################################################################
# HISTORY:
# 2011-03-27
# o Created.
############################################################################
