#########################################################################/**
# @set "class=RGData"
# @RdocMethod fitIWPCA
#
# @title "Robust fit of a line through a multidimensional data"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{N-times-K @matrix where N is the number of observations and
#           K is the number of dimensions.}
#  \item{w}{Weights for all data points.}
#  \item{constraint}{@character specifying which additional contraint to 
#   be used to specify the offset parameters along the fitted line.
#   If \code{"diagonal"}, the offset vector will be point on the line 
#   that is closest to the diagonal line (1,...,1).
#   If \code{"max"}, the offset vector will the point on the line that is
#   as "great" as possible, but still such that each of its components is
#   less than the corresponding minimal signal. This will guarantee that
#   no negative signals are created in the backward transformation.}
#  \item{...}{Additional arguments accepted by @see "aroma.light::iwpca".
#   For instance, a N @vector of weights for each (observation) may be
#    given, otherwise they get the same weight.}
#  \item{verbose}{If @TRUE, detailed information is printed during fit.}
# }
#
# \value{
#   Returns a @list contain the fitted parameter etc. [MORE].
# }
#
# \details{
#   This method uses re-weighted principal component analysis (IWPCA)
#   to fit a the nodel \eqn{y_n = a + bx_n + eps_n} where \eqn{y_n},
#   \eqn{a}, \eqn{b}, and \eqn{eps_n} are vector of the K and \eqn{x_n}
#   is a scalar. 
#
#   The algorithm is:
#    For iteration i
#    1) Fit a line \eqn{L} through the data close using weighted PCA
#       with weights \eqn{\{w_n\}}. Let
#         \eqn{r_n = \{r_{n,1},...,r_{n,K}\}}
#       be the \eqn{K} principal components.
#    2) Update the weights as
#         \eqn{w_n <- 1 / \sum_{2}^{K} (r_{n,k} + \epsilon_r)}
#       where we have used the residuals of all but the first principal
#       component.
#    3) Find the point a on \eqn{L} that is closest to the
#       line \eqn{D=(1,1,...,1)}. Similarily, denote the point on D that is 
#       closest to \eqn{L} by \eqn{t=a*(1,1,...,1)}.
# }
#
# @author
#
# %examples "RGData.fitMultiIWPCA.Rex"
#
# \seealso{
#   This is an internal method used by the @seemethod "fitMultiscanAffine"
#   method, which in addiion uses @seemethod "fitPairIWPCA".
#   Internally the function @see "aroma.light::iwpca" is used to fit a line 
#   through the data cloud and 
#   the function @see "aroma.light::distanceBetweenLines" to find the closest
#   point to the diagonal (1,1,...,1).
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("fitIWPCA", "RGData", function(static, X, w=NULL, constraint=c("diagonal", "max"), ..., aShift=rep(0, ncol(X)), Xmin=NULL, bootstrap=FALSE, R=6, verbose=FALSE) {
  # A function which when applied to data returns a vector containing the
  # statistic(s) of interest. 
  # This function should be such that it can be passed to 
  # boot(..., statistic=statistic) in the 'boot' package. For this 
  # statistic() must take at least two arguments. The first argument passed
  # will always be the original data. The second will be a vector of 
  # indices.
  bootCount <- -1;  # First one will always be the original data.
  statistic <- function(X, w=NULL, index=seq(nrow(X)), ..., constraint="diagonal", Xmin=NULL, aShift=rep(0, ncol(X))) {
    gc();
    bootCount <<- bootCount + 1;
    if (verbose == TRUE && bootCount > 0) 
      cat(sprintf("Bootstrap %d:\n", as.integer(bootCount)));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Resample the orginal data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    X <- X[index,];
    if (!is.null(w))
      w <- w[index];

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit an K-dimensional line through the data using iterative
    # re-weighted PCA.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fit <- iwpca(X, w=w, ...);                          # {aroma.light}

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get the fitted line L
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get the center of the fitted line...
    ax <- fit$xMean;
    names(ax) <- NULL;
  
    # ...and the fitted eigenvectors (u1,u2,...,uK) 
    # with ui*uj = 0; i!=j and ui*ui = 1.
    U <- t(fit$vt);
    colnames(U) <- rownames(U) <- NULL;
  
    # The fitted scale parameters b=(b[1],b[2],...,b[K]) where
    # the elements are rescaled such that b[1] == 1.
    # [ min(b[i]) == 1. Before it was such that b[1] == 1, but this
    # is probably better. (Indeed not; this is not good if one do more
    # than one estimate per array, e.g. printtip etc. /HB 2004-04-26) ]
    U1 <- U[,1];
    bx <- as.vector(U1/U1[1]);

    # Shift the data. 
    # [ This is for instance useful if fitting towards the diagonal line
    #   and resampling under H0: y_i = alpha + z_i and /040102 ]
    ax <- ax + aShift;

    if (identical(constraint, "diagonal")) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Find the point t on the fitted line that is closest to the
      # points s on the "diagonal" line (1,1,...,1) in K-space.
      # [This works also for lines in two dimension.]
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # x(s) is the fitted line (the first IWPCA component)
      # y(t) is the diagonal line
      ay <- rep(0,length(ax));      # (0,0,...,0)
      by <- rep(1,length(ay));      # (1,1,...,1)
  
      dbl <- distanceBetweenLines(ax=ax,bx=bx, ay=ay,by=by);  # {aroma.light}
      a     <- as.vector(dbl$xs);
      adiag <- as.vector(dbl$yt);
    } else if (identical(constraint, "max")) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Find the "greatest" point t on the fitted line that is within
      # the cube C whose upper limits are defined by the minimum value
      # in each channel.
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Find the minimal value of each X component.
      if (is.null(Xmin))
        Xmin <- apply(X, MARGIN=2, FUN=min, na.rm=TRUE);
      # For each component k, find the value t such that 
      #    ax[k] + bx[k]*t[k] == Xmin[k]  <=> t[k] == (Xmin[k] - ax[k])/bx[k]
      t <- (Xmin-ax)/bx;
      # Choose minimum t[k]
      # Now, amax <<= Xmin if amax <- ax - bx[k]*min(t) where <<= (\prec) 
      # means componentswise less or equal than.
      a <- ax + bx*min(t);
      adiag <- rep(0, length(ax));
    }

    # Return the statistic
    t <- c(a=a);
    t <- c(t, b=bx);
    t <- c(t, adiag=adiag);
    t <- c(t, U=as.vector(U));    
    t <- c(t, niter=fit$nbrOfIterations * (fit$converged*2-1));
    t;
  } # statistic()

  if (!is.matrix(X))
    throw("Argument 'X' must be a matrix.");
  N <- nrow(X);
  K <- ncol(X);
  if (K == 1)
    throw("Argument 'X' must have two or more columns: ", K);
  if (N < K)
    throw("Argument 'X' must have at least as many rows as columns: ", N, " < ", K);

  constraint <- match.arg(constraint);

  if (is.null(aShift)) {
    aShift <- rep(0, ncol(X));
  }

  if (identical(constraint, "max")) {
     Xmin <- apply(X, MARGIN=2, FUN=min, na.rm=TRUE);
  } else {
     Xmin <- NULL;
  }

  require("R.basic") || throw("Package 'R.basic' not found.");

  if (bootstrap == TRUE)  {
    require("boot") || 
        throw("Package 'boot' not found. Can perform bootstrap resampling.");
    if (is.null(R))
      R <- 50;
    boot <- boot(X, statistic=statistic, w=w, R=R, constraint=constraint, Xmin=Xmin, aShift=aShift, ...);
    # Save memory by deleting obsolete large elements
    t0 <- boot$t0;
    t  <- boot$t;
    colnames(t) <- names(t0);
    rm(boot);
  } else {
    nresamples <- 1;
    t0 <- statistic(X, w=w, index=seq(nrow(X)), constraint=constraint, Xmin=Xmin, aShift=aShift, ...);
    t  <- NULL;
  }

  a <- t0[regexpr("^a[0-9]*$", names(t0)) != -1];
  b <- t0[regexpr("^b[0-9]*$", names(t0)) != -1];
  adiag <- t0[regexpr("^adiag[0-9]*$", names(t0)) != -1];
  U <- t0[regexpr("^U[0-9]*$", names(t0)) != -1];
  U <- matrix(U, nrow=sqrt(length(U)));
  niter <- as.integer(abs(t0["niter"]));
  converged <- (niter > 0);

  AffineModelFit(t0=t0, t=t, a=a, b=b, adiag=adiag, eigen=U, y=X, converged=converged, nbrOfIterations=niter);
}, static=TRUE, protected=TRUE, deprecated=TRUE); # fitIWPCA()




#########################################################################/**
# @RdocMethod fitMultiscanAffine
#
# @title "Fits an affine model to signals from a multi-scanned slide"
#
# \description{
#   @get "title".
#
#   This function is used internally by @seemethod "calibrateMultiscan" and
#   is not intended to be used by the end-user. Developers of new algorithms
#   may use it, but there is no guarantee that it will exist in the future.
# }
#
# @synopsis
#
# \arguments{
#  \item{slides}{Slides to be used in the fit \emph{and} that are to be 
#   calibrated. If @NULL, all slides are considered.}
#  \item{channels}{Specifies which channels to be fitted. If @NULL, all
#   channels are considered.}
#  \item{method}{@character specifying method for IWPCA fitting. 
#   If \code{"L1"}, the distances are minimised in \eqn{L_1}. 
#   No other methods are available.}
#  \item{maxIter}{The maximum number of re-weigthed iterations.}
#  \item{acc}{The accuracy for the stopping criterion.}
#  \item{satSignal}{Signals equal to or above this threshold will not
#    be used in the fitting.}
#  \item{bootstrap}{If @TRUE, bootstrap is applied to estimate the
#    mean and standard deviation of the bootstraped estimates.}
#  \item{R}{The number of bootstrap samples. Only effective if 
#   \code{bootstrap==TRUE}.}
#  \item{reps}{Small value, which is added to the absolute value of 
#    the residuals to avoid "1/0" weights, i.e. "(1/0+reps)".}
#  \item{verbose}{If @TRUE, extra information is printed while running.}
#  \item{...}{Additional arguments accepted by @see "aroma.light::iwpca".
#   For instance, a N @vector of weights for each observation may be
#    given, otherwise they get the same weight.}
# }
#
# \value{
#   Returns a @list of fitted @see "AffineModelFit" objects for
#   each channel.
# }
#
# \details{
#  Fitting is done by iterated re-weighted principal component analysis
#  (IWPCA).
# }
#
# @author
#
# \references{
#   [1] @include "../incl/BengtssonH_etal_2004.bib.Rdoc" \cr
# }
#
# \seealso{
#   To calibrate the data using these fits see @seemethod "calibrateMultiscan".
#   For a similar strategy for between channel normalization see
#   @seemethod "normalizeAffine".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("fitMultiscanAffine", "RGData", function(this, slides=NULL, channels=c("R", "G"), aShift=NULL, method="L1", maxIter=30, reps=0.02, acc=1e-4, satSignal=2^16-1, bootstrap=FALSE, R=NULL, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  slides <- validateArgumentSlides(this, slides=slides);
  method <- match.arg(method);

  if (reps < 0)
    throw("Argument 'reps' must be non-negative: ", reps);

  # Number of scans
  nscans <- length(slides);
  if (nscans == 1) {
    throw("Can not fit affine multiscan model. Argument 'slides' must contain at least two slides: ", nscans);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Fit the model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fits <- list();

  # Fit the model for each channel seperately
  for (channel in channels) {
    X <- this[[channel]][,slides];

    # Use only finite and non-saturated observations.
    ok <- apply(X, MARGIN=1, FUN=function(x) all(is.finite(x)));
    X <- X[ok,];
    ok <- apply(X, MARGIN=1, FUN=function(x) all(x < satSignal));
    X <- X[ok,];
    rm(ok);

    # Number of observations
    nobs <- nrow(X);

    if (verbose) 
      cat(sprintf("Estimating bias based on %d observations.\n", nobs));

    fit <- RGData$fitIWPCA(X, aShift=aShift, maxIter=maxIter, acc=acc, reps=reps, bootstrap=bootstrap, R=R, ...);

    # X os not needed anymore, because it is a subset of all spots!
    rm(X);

    fit$slides <- slides;
    fits[[channel]] <- fit;
  } # for (channel ...)

  fits;
}, protected=TRUE) # fitMultiscanAffine()


############################################################################
# HISTORY:
# 2005-02-02
# o Excluded fitIWPCA() and fitMultiscanAffine() from RGData() since they
#   are now in the matrix class.
############################################################################
