setConstructorS3("RawCopyNumbers", function(cn=NULL, x=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cn':
  if (!is.null(cn)) {
    if (!is.vector(cn)) {
      throw("Argument 'cn' must be a vector: ", mode(cn)[1]);
    }

    if (!is.numeric(cn)) {
      throw("Argument 'cn' must be a numeric: ", class(cn)[1]);
    }
  }

  if (!is.null(x)) {
    if (!is.vector(x)) {
      throw("Argument 'x' must be a vector: ", mode(x)[1]);
    }

    if (!is.numeric(x)) {
      throw("Argument 'x' must be a numeric: ", class(x)[1]);
    }

    n <- length(cn);
    if (length(x) != n) {
      throw("Argument 'x' and 'cn' are of different lengths: ", length(x), " != ", n);
    }
  }

  extend(Object(), "RawCopyNumbers", 
    cn = cn,
    x = x
  )
})


setMethodS3("as.character", "RawCopyNumbers", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Number of loci: %d", nbrOfLoci(this)));
  fields <- getLociFields(this);
  fields <- sapply(fields, FUN=function(field) {
    values <- this[[field]];
    mode <- mode(values);
    sprintf("%s [%dx%s]", field, length(values), mode);
  })
  fields <- paste(fields, collapse=", ");
  s <- c(s, sprintf("Loci fields: %s", fields));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE) 


setMethodS3("nbrOfLoci", "RawCopyNumbers", function(this, ...) {
  length(this$x);
})

setMethodS3("getPhysicalPositions", "RawCopyNumbers", function(this, ...) {
  this$x;
})

setMethodS3("getCNs", "RawCopyNumbers", function(this, ...) {
  this$cn;
})

setMethodS3("as.data.frame", "RawCopyNumbers", function(x, ...) {
  # To please R CMD check
  this <- x;

  data.frame(x=this$x, cn=this$cn);
})

setMethodS3("summary", "RawCopyNumbers", function(object, ...) {
  # To please R CMD check
  this <- object;

  summary(as.data.frame(this));
})


setMethodS3("getLociFields", "RawCopyNumbers", function(this, ...) {
  c("cn", "x");
})

setMethodS3("sort", "RawCopyNumbers", function(x, ...) {
  # To please R CMD check
  this <- x;

  res <- clone(this);
  o <- order(res$x);
  for (field in getLociFields(res)) {
    res[[field]] <- res[[field]][o];
  }
  res;
})

setMethodS3("getXY", "RawCopyNumbers", function(this, sort=TRUE, ...) {
  xy <- data.frame(x=this$x, y=this$cnl);
  if (sort)
    xy <- xy[order(xy$x),];
  xy;
})

setMethodS3("plot", "RawCopyNumbers", function(x, xlab="Physical position", ylab="Relative copy number", ylim=c(-3,3), pch=20, ...) {
  # To please R CMD check
  this <- x;

  plot(this$x, this$cn, ylim=ylim, xlab=xlab, ylab=ylab, pch=pch, ...);
})


setMethodS3("points", "RawCopyNumbers", function(x, pch=20, ...) {
  # To please R CMD check
  this <- x;

  points(this$x, this$cn, pch=pch, ...);
})

setMethodS3("lines", "RawCopyNumbers", function(x, ...) {
  # To please R CMD check
  this <- x;

  x <- this$x;
  o <- order(x);
  x <- x[o];
  y <- this$cn[o];
  lines(x, y, ...);
})


setMethodS3("gaussianSmoothing", "RawCopyNumbers", function(this, sd=10e3, ...) {
  # Get a sorted copy
  res <- sort(this);

  res$cn <- gaussianSmoothing(res$cn, x=res$x, sd=sd, na.rm=TRUE, ...);

  res;
})


setMethodS3("extractRawCNs", "default", function(...) {
  extractRawCopyNumbers(...);
})



setMethodS3("extractRawCopyNumbers", "profileCGH", function(object, ...) {
  pv <- object$profileValues;
  RawCopyNumbers(cn=pv$LogRatio, x=pv$PosBase);
})


setMethodS3("extractRawCopyNumbers", "DNAcopy", function(object, ...) {
  data <- object$data;
  RawCopyNumbers(cn=data[[3]], x=data$maploc);
})




############################################################################
# HISTORY:
# 2007-08-22
# o Created.  Need a generic container for holding copy number data and
#   to plot them nicely.
############################################################################
