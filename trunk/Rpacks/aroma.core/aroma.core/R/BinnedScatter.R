setMethodS3("binScatter", "matrix", function(x, nbin=128, orderBy="density", decreasing=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'x':
  dim <- dim(x);
  if (dim[2] != 2) {
    throw("Argument 'x' must be a two-column matrix: ", dim[2]);
  }

  # Argument 'orderBy':
  if (!is.null(orderBy)) {
    orderBy <- match.arg(orderBy, c("density"));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate the (x,y) density
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ok <- whichVector(is.finite(x[,1]) & is.finite(x[,2]));
  x <- x[ok,,drop=FALSE];
  rm(ok);
  map <- geneplotter:::.smoothScatterCalcDensity(x, nbin=nbin);

  xm <- map$x1;
  ym <- map$x2;
  dens <- map$fhat;

  nx <- length(xm);
  ny <- length(ym);
  dx <- x[,1]-xm[1];
  dy <- x[,2]-ym[1];
  w <- xm[nx] - xm[1];
  h <- ym[ny] - ym[1];
  ixm <- round(dx/w * (nx - 1));
  iym <- round(dy/h * (ny - 1));
  idens <- dens[1 + iym*nx + ixm];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup return structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list(
    data=x,
    density=idens,
    map=map,
    params=list(nbin=nbin, orderBy=orderBy, decreasing=decreasing)
  );

  class(res) <- c("BinnedScatter", class(res));


  # Reorder data?
  if (!is.null(orderBy)) {
    o <- order(res[[orderBy]], decreasing=decreasing);
    res$data <- res$data[o,,drop=FALSE];
    res$density <- res$density[o];
  }

  res;
}) # calculateScatterDensity()


setMethodS3("points", "BinnedScatter", function(object, ...) {
  points(object$data, ...);
})


setMethodS3("plot", "BinnedScatter", function(object, ...) {
  plot(object$data, ...);
})


setMethodS3("subset", "BinnedScatter", function(object, subset, ...) {
  object$data <- object$data[subset,,drop=FALSE];
  object$density <- object$density[subset];
  object;
})


setMethodS3("subsample", "BinnedScatter", function(object, size=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  weightFcn <- function(object, ...) {
    w <- 1/object$density;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  n <- dim(object$data)[1];

  # Argument 'size':
  if (is.null(size)) {
    size <- n;
  } else {
    size <- Arguments$getNumeric(size, range=c(0, n));
    if (size < 1) {
      size <- round(size*n);
      if (size > n)
        size <- n;
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate sample weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  w <- weightFcn(object);

  # Standarize weights
  w <- w / sum(w, na.rm=TRUE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Randomized sampling according to weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subset <- sample(1:length(w), size=size, prob=w, replace=FALSE);
  rm(w);
  res <- subset(object, subset=subset, ...);

  res;
}) # subsample()



############################################################################
# HISTORY:
# 2008-11-14
# o Created.
############################################################################ 
