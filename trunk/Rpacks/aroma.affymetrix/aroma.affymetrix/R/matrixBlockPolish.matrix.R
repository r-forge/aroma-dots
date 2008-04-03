###########################################################################/**
# @set "class=matrix"
# @RdocMethod matrixBlockPolish
#
# @title "Applies a polishing function to blocks of rows and columns repeatedly"
#
# \description{
#  @get "title".
# }
# 
# @synopsis
#
# \arguments{
#   \item{y}{A @numeric KxN @matrix.}
#   \item{blockSizes}{A positive @integer @vector of length two.}
#   \item{FUN}{A @function taking @numeric arguments \code{y} and 
#      \code{x} and returns a @numeric object with either a scalar
#      or the same number of elements as in \code{y}.}
#   \item{...}{Additional arguments passed to the \code{FUN} @function.}
#   \item{tol}{A positive threshold specifying when the algorithm has
#      converged.}
#   \item{maxIter}{The maximum number of iterations.}
#   \item{returnEffects}{If @TRUE, the row and column effects are returned,
#      otherwise not.}
# }
#
# \value{
#   Returns a named @list.
# }
#
# @examples "../incl/matrixBlockPolish.matrix.Rex"
#
# @author
#
# \seealso{
#  @see "stats::medpolish".
#  @see "aroma.light::medpolish".
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("matrixBlockPolish", "matrix", function(y, blockSizes=c(1,1), FUN, ..., tol=0.01, maxIter=10, returnEffects=FALSE) {
  # Argument 'blockSizes':
  blockSizes <- rep(as.integer(blockSizes), length.out=2);

  # Argument 'FUN':
  if (!is.function(FUN)) {
    throw("Argument 'FUN' is not a function: ", class(FUN)[1]);
  }

  ranges <- vector("list", length=2);
  for (dd in 1:2) {
    idxs1 <- seq(from=1, to=dim(y)[dd], by=blockSizes[dd]);
    idxs2 <- c(idxs1[-1]-1, dim(y)[dd]);
    ranges[[dd]] <- cbind(from=idxs1, to=idxs2);
  }
  rm(idxs1, idxs2);

  if (returnEffects) {
    blockSizes <- sapply(ranges, FUN=function(r) r[,2]-r[,1]+1);
    maxBlockSizes <- sapply(blockSizes, FUN=max);
    effects <- vector("list", length=2);
    names(effects) <- c("rows", "columns");
    for (dd in 1:2) {
      nbrOfBlocks <- nrow(ranges[[dd]]);
      effects[[dd]] <- matrix(as.double(NA), nrow=nbrOfBlocks, ncol=maxBlockSizes[dd]*dim(y)[3-dd]);
    }
  }

  x <- matrix(seq(length=length(y)), nrow=nrow(y), ncol=ncol(y));
  oldSum <- 0;
  for (ii in seq(length=maxIter)) {
    for (dd in 1:2) {
      range <- ranges[[dd]]
      froms <- range[,1]; 
      tos <- range[,2];
      nbrOfBlocks <- length(froms);
  
      for (kk in seq(length=nbrOfBlocks)) {
        idxs <- froms[kk]:tos[kk];
    
        # Get data
        if (dd == 1) {
          xB <- x[idxs,,drop=FALSE];
          yB <- y[idxs,,drop=FALSE];
        } else if (dd == 2) {
          xB <- x[,idxs,drop=FALSE];
          yB <- y[,idxs,drop=FALSE];
        }
    
        # Polish data
        yB2 <- FUN(yB, xB, ...);
        rm(xB);
    
        yB <- yB - yB2;
  
        if (returnEffects) {
          effects[[dd]][kk,] <- yB2;
        }

        rm(yB2);
    
        # Update data
        if (dd == 1) {
          y[idxs,] <- yB;
        } else if (dd == 2) {
          y[,idxs] <- yB;
        }
    
        rm(yB);
      } # for (kk ...)
    } # for (dd ...)

    newSum <- sum(abs(y), na.rm=TRUE);
    converged <- (newSum == 0 || abs(newSum - oldSum) < tol * newSum);
    if (converged)
      break;

    oldSum <- newSum;
  } # for (ii ...)

  res <- list(residuals=y, converged=converged, iter=ii);
  if (returnEffects) {
    effects <- lapply(effects, FUN=drop);
    res <- c(list(row=effects[[1]], col=effects[[2]]), res);
  }
  class(res) <- c("matrixBlockPolish");

  res;
})



############################################################################
# HISTORY: 
# 2008-04-02
# o Verified against median polish.
# o Created.
############################################################################ 
