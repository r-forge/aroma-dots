setMethodS3("lapplyByGene", "RGArray", function(object, FUN, ...) {
  X <- object;

  names <- rownames(X);
  groupBy <- by(seq(along=names), INDICES=names, FUN=function(x) x);

  res <- base::lapply(groupBy, FUN=function(rows) {
    Xgg <- X[rows,,,drop=FALSE];
    FUN(Xgg, ...);
  });

  res;
}, protected=TRUE)

#
# \example{\dontrun{
#  Xm <- averageByGene(X);
#  Xm2 <- averageByGene(Xm);
#  # Sanity check
#  stopifnot(identical(Xm2, Xm));
# }}
#
setMethodS3("averageByGene", "RGArray", function(object, FUN=mean, na.rm=TRUE, ...) {
  # Argument 'FUN':
  stopifnot(is.function(FUN));

  fcn <- function(Xg, na.rm=TRUE, ...) {
    apply(Xg, MARGIN=c(2,3), FUN=FUN, na.rm=na.rm, ...);
  } # fcn()

  # Get a list of length J consisting of CxI matrices
  XavgList <- lapplyByGene(object, FUN=fcn, na.rm=na.rm, ...);

  # Turn into a CxIxJ array
  Xavg <- abind::abind(XavgList, along=3L);

  # Permute to a JxCxI array
  Xavg <- aperm(Xavg, perm=c(3,1,2));

  # Make into an RGArray (of the same class)
  class(Xavg) <- class(object);

  Xavg;
}) # averageByGene()


setMethodS3("swapRG", "RGArray", function(this, ...) {
  cc <- seq(length=dim(this)[2]);
  this[,cc,] <- this[,rev(cc),];
  this;
}) # swapRG()


##############################################################################
# HISTORY:
# 2011-12-06
# o Added swapRG() for RGArray.
# 2011-09-14
# o Added averageByGene() for RGArray.  Requires package 'abind'.
# o Added lapplyByGene() for RGArray.
# o Created.
##############################################################################
