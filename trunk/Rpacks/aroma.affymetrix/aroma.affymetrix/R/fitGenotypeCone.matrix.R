###########################################################################/**
# @set "class=matrix"
# @RdocMethod fitGenotypeCone
#
# @title "Fits an affine transformation to allele A and allele B data"
#
# \description{
#  @get "title" using robust estimators.
# }
# 
# @synopsis
#
# \arguments{
#   \item{y}{A @numeric Nx2 @matrix with one column for each allele and
#      where N is the number of data points.}
#   \item{alpha}{A @numeric @vector of decreasing values in (0,1).
#      This parameter "determines how far we are willing to press the
#      boundary of the [genotype cone]".  Lowering \code{alpha} expand
#      the cone.  When \code{alpha} goes to zero, all data points will
#      be on or inside the cone.}
#   \item{q,Q}{Percentiles in [0,100] for which data points that are 
#      below (above) will be assigned zero weight in the fitting of
#      the parameters.}
#   \item{...}{Additional arguments passed to @see "sfit::cfit".}
# }
#
# \value{
#   Returns the parameter estimates as a named @list with elements:
#   \itemize{
#    \item{M}{An estimate of the three vertices defining the genotype
#      triangle.  These three vertices are describes as an 2x3 @matrix
#      with column \code{origin}, \code{AA}, and \code{BB}.}
#    \item{Minv}{The inverse of \code{M}.}
#    \item{origin}{The estimate of the shift.}
#    \item{W}{The estimate of shear/rotation matrix with columns 
#             \code{AA} and \code{BB}.}
#    \item{Winv}{The inverse of \code{W}.}
#    \item{params}{The parameters used for the fit, i.e. 
#       \code{alpha}, \code{q}, \code{Q}, and  those passed in \code{...}.}
#    \item{dimData}{The dimension of the input data.}
#   }
# }
#
# @examples "../incl/fitGenotypeCone.matrix.Rex"
#
# @author
#
# \seealso{
#  To backtransform data fitted using this method, 
#  see @see "backtransformGenotypeCone".
#  Internally @see "sfit::cfit" of the \pkg{sfit} package is used.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("fitGenotypeCone", "matrix", function(y, alpha=c(0.10, 0.075, 0.05, 0.03, 0.01), q=2, Q=98, ...) {
  require("sfit") || throw("Package 'sfit' not found.");

  # Fit simplex of (y_A,y_B)
  fit <- sfit::cfit(y, alpha=alpha, q=q, Q=Q, ...);

  M <- fit$M;
  colnames(M) <- c("A", "B");
  clazz <- class(M);

  # Re-arrange vertices in the order (origin, AA, BB)
  idxOrigin <- which.min(apply(M, MARGIN=1, FUN=function(u) sum(u^2)));
  origin <- M[idxOrigin,];
  M <- M[-idxOrigin,];
  idxBBAA <- order(apply(M, MARGIN=1, FUN=function(u) diff(u)));
  M <- M[idxBBAA,];
  M <- rbind(origin, M);
  rownames(M) <- c("origin", "AA", "BB");
  class(M) <- clazz;
  fit$M <- M;

  W <- M[c("AA","BB"),];
  W <- t(W) - origin;
  W <- W / W[1,1];
 
  # Find the inverse
  Winv <- solve(W);

  Minv <- t(Winv %*% (t(M)-origin));
  class(Minv) <- clazz;
  fit$Minv <- Minv;

  W <- t(W);
  Winv <- t(Winv);

  fit$origin <- origin;
  fit$W <- W;
  fit$Winv <- Winv;

##  Excluded. /HB 2007-06-12
##  # Re-arrange X too
##  if (!is.null(fit$X)) {
##    fit$X <- fit$X[,oB];  ### 'oB'??? HB 2007-06-11
##  }

  fit$params <- list(
    alpha=alpha,
    q=q,
    Q=Q,
    ...
  );

  fit$dimData <- dim(y);

  fit;
}, private=TRUE) # fitGenotypeCone()



############################################################################
# HISTORY:
# 2008-02-14
# o Added a self-contained example for fitGenotypeCone().
# 2007-09-08
# o Added 'dimData' to the return structure.
# 2007-06-12
# o Commented the code for re-arranging fit$X (only if retX=TRUE).
#   Code not really needed since backtransformGenotypeCone() is used.
# 2007-06-11
# o Calls sfit::cfit() explicitly.
# 2007-06-04
# o Added Rdoc comments.
# 2006-05-08
# o Created.
############################################################################
