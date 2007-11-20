###########################################################################/**
# @RdocDefault normalizeFragmentLength
#
# @title "Normalizes signals for PCR fragment-length effects"
#
# \description{
#  @get "title". Some or all signals are used to estimated the 
#  normalization function.  All signals are normalized.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{y}{A @numeric @vector of length K of signals to be normalized
#     across E enzymes.}
#   \item{fragmentLengths}{An @integer KxE @matrix of fragment lengths.}
#   \item{targetFcns}{A @list of E @functions - one per enzyme.}
#   \item{subsetToFit}{The subset of data points used to fit the 
#      normalization function.
#      If @NULL, all data points are considered.}
#   \item{.isLogged}{A @logical.}
#   \item{...}{Additional arguments passed to @see "stats::lowess".}
#   \item{.returnFit}{A @logical.}
# }
#
# \value{
#   Returns a @numeric @vector of the normalized signals.
# }
#
# \section{Multi-enzyme normalization}{
#  It is assumed that the fragment-length effects from multiple enzymes
#  added (with equal weights) on the intensity scale.
#  The fragment-length effects are fitted for each enzyme separately based
#  on units that are exclusively for that enzyme. 
#  \emph{If there are no or very such units for an enzyme, the assumptions 
#  of the model are not met and the fit will fail with an error.}
#  Then, from the above single-enzyme fits the average effect across
#  enzymes is the calculated for each unit that is on multiple enzymes.
# }
#
# \examples{
#   @include "../incl/normalizeFragmentLength-ex1.Rex"
#
#   @include "../incl/normalizeFragmentLength-ex2.Rex"
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("normalizeFragmentLength", "default", function(y, fragmentLengths, targetFcns=NULL, subsetToFit=NULL, .isLogged=TRUE, ..., .returnFit=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  y <- Arguments$getDoubles(y, disallow=NULL);
  
  # Argument 'fragmentLengths':
  if (!is.matrix(fragmentLengths)) {
    if (is.vector(fragmentLengths)) {
      fragmentLengths <- as.matrix(fragmentLengths);
    } else {
      throw("Argument 'fragmentLengths' must be a matrix: ", 
                                                class(fragmentLengths)[[1]]);
    }
  }
  nbrOfEnzymes <- ncol(fragmentLengths);
  allEnzymes <- seq(length=nbrOfEnzymes);
  for (ee in allEnzymes) {
    fragmentLengths[,ee] <- Arguments$getDoubles(fragmentLengths[,ee], length=length(y), disallow=NULL);
  }

  # Argument 'targetFcns':
  if (!is.null(targetFcns)) {
    if (!is.list(targetFcns)) {
      if (nbrOfEnzymes == 1) {
        targetFcns <- list(targetFcns);
      } else {
        throw("Argument 'targetFcns' is not a list: ", class(targetFcns)[1]);
      }
    }
    if (length(tagetFcns) != nbrOfEnzymes) {
      throw("Number of elements in 'targetFcns' does not match the number of columns in 'fragmentLengths': ", length(tagetFcns), " != ", nbrOfEnzymes);
    }

    # Validate each element
    for (ee in allEnzymes) {
      if (!is.function(targetFcns[[ee]])) {
        throw("One element in 'targetFcns' is not a function: ", class(targetFcns[[ee]])[1]);
      }
    }
  }
  
  # Argument 'subsetToFit':
  if (!is.null(subsetToFit)) {
    subsetToFit <- Arguments$getIndices(subsetToFit, length=length(y));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate normalization function and predict the signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit smooth curve to each enzyme separately
  hasFL <- is.finite(fragmentLengths);
  # Count the number of enzymes per units
  countFL <- rep(0, nrow(hasFL));
  for (ee in allEnzymes)
    countFL <- countFL + as.integer(hasFL[,ee]);
  isSingleEnzymed <- (countFL == 1);

  okY <- is.finite(y);

  # KxE matrix of corrections
  dy <- matrix(NA, nrow=nrow(fragmentLengths), ncol=nbrOfEnzymes);

  if (.returnFit) {
    fits <- vector("list", nbrOfEnzymes);
  }

  for (ee in allEnzymes) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (a) Fit normalization function
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit only to units with known length and non-missing data points.
    ok <- (hasFL[,ee] & isSingleEnzymed & okY);

    # Sanity check
    if (sum(ok) == 0) {
      throw("Cannot fit normalization function to enzyme, because there are no (finite) data points that are unique to this enzyme: ", ee);
    }

    if (!is.null(subsetToFit)) {
      ok[-subsetToFit] <- FALSE;
    }

    # Sanity check
    if (sum(ok) == 0) {
      throw("Cannot fit normalization function to enzyme, because there are no (finite) data points that are unique to this enzyme for the subset requested: ", ee);
    }

    fl <- fragmentLengths[,ee];
    suppressWarnings({
      fit <- lowess(fl[ok], y[ok], ...);
      class(fit) <- "lowess";
    })
    rm(ok);

    if (.returnFit)
      fits[[ee]] <- fit;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (b) Calculate correction factor
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dy[,ee] <- predict(fit, newdata=fl);

    # Normalize toward a target function?
    if (!is.null(targetFcns)) {
      yTargetPred <- targetFcns[[ee]](fl);
      dy[,ee] <- dy[,ee] - yTargetPred;
      rm(yTargetPred);
    }

    if (.returnFit)
      fits[[ee]] <- list(fit=fit, dy=dy);

    rm(fit, fl);

  } # for (ee ...)
  rm(hasFL, isSingleEnzymed);

  # Sum on the non-log scale.
  if (.isLogged) {
    dy <- 2^dy;
    dy <- rowSums(dy, na.rm=TRUE);
    dy <- dy / countFL;
    dy <- log2(dy);
  }
  rm(countFL);

  # Transform signals
  y[okY] <- y[okY] - dy[okY];
  

  if (.returnFit)
    attr(y, "modelFit") <- fits;

  y;
}, private=TRUE)


############################################################################
# HISTORY:
# 2007-11-19
# o Added Rdoc examples. From these simulation examples, it looks like the
#   multi-enzyme normalization method works.
# o Updated normalizeFragmentLength() to handle multiple enzymes.
# 2006-11-28
# o Created.
############################################################################
