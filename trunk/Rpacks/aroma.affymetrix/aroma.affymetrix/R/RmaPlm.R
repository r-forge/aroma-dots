###########################################################################/**
# @RdocClass RmaPlm
#
# @title "The RmaPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the log-additive model used in RMA.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeLevelModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#   For a single unit group, the log-additive model of RMA is:
#
#    \deqn{log_2(y_{ij}) = \beta_i + \alpha_j + \varepsilon_{ij}}
#
#   where \eqn{\beta_i} are the chip effects for arrays \eqn{i=1,...,I}, 
#   and \eqn{\alpha_j} are the probe affinities for probes \eqn{j=1,...,J}.
#   The \eqn{\varepsilon_{ij}} are zero-mean noise with equal variance.
#   The model is constrained such that \eqn{\sum_j{\alpha_j} = 0}.
#
#   Note that all PLM classes must return parameters on the intensity scale.
#   For this class that means that \eqn{\theta_i = 2^\beta_i} and 
#   \eqn{\phi_i = 2^\alpha_i} are returned.
# }
#
# @author
#*/###########################################################################
setConstructorS3("RmaPlm", function(..., tags="RMA") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(affyPLM) || throw("Package 'affyPLM' not loaded.");

  extend(ProbeLevelModel(..., tags=tags), "RmaPlm")
})



setMethodS3("getProbeAffinities", "RmaPlm", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the probe affinities (and create files etc)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  paf <- NextMethod("getProbeAffinities", this, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the encode and decode functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setEncodeFunction(paf, function(groupData, ...) {
    phi <- .subset2(groupData, "phi");
    stdvs <- .subset2(groupData, "sdPhi");
    outliers <- .subset2(groupData, "phiOutliers");

    # Encode outliers as the sign of 'pixels'; -1 = TRUE, +1 = FALSE
    pixels <- sign(0.5 - as.integer(outliers));

    list(intensities=phi, stdvs=stdvs, pixels=pixels);
  })

  setDecodeFunction(paf,  function(groupData, ...) {
    intensities <- .subset2(groupData, "intensities");
    stdvs <- .subset2(groupData, "stdvs");
    pixels <- .subset2(groupData, "pixels");

    # Outliers are encoded by the sign of 'pixels'.
    outliers <- as.logical(1-sign(pixels));

    list(
      phi=intensities, 
      sdPhi=stdvs, 
      phiOutliers=outliers
    );
  })

  paf;
})
  


###########################################################################/**
# @RdocMethod getFitFunction
#
# @title "Gets the low-level function that fits the PLM"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @function.
# }
#
# \author{
#   Henrik Bengtsson and Ken Simpson (WEHI) utilizing Ben Bolstad's 
#   \pkg{affyPLM} package.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitFunction", "RmaPlm", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # rmaModel02()
  # Author: Henrik Bengtsson, UC Berkeley. 
  # Requires: affyPLM() by Ben Bolstad.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rmaModel01 <- function(y, psiCode=0, psiK=1.345){
    # Assert right dimensions of 'y'.
    if (length(dim(y)) != 2) {
      str(y);
      stop("Argument 'y' must have two dimensions: ", 
                                                paste(dim(y), collapse="x"));
    }

    # Log-additive model
    y <- log(y, base=2)

    # Fit model using affyPLM code
    fit <- .Call("R_rlm_rma_default_model", y, psiCode, psiK, PACKAGE="affyPLM");

    # Extract probe affinities and chip estimates
    J <- ncol(y);
    I <- nrow(y);
    est <- fit$Estimates;
    se <- fit$StdErrors;

    # Chip effects
    beta <- est[1:J];

    # Probe affinities
    alpha <- est[(J+1):length(est)];
    alpha[length(alpha)] <- -sum(alpha[1:(length(alpha)-1)]);

    # Estimates on the intensity scale
    theta <- 2^beta;
    phi <- 2^alpha;

    # The RMA model is fitted with constraint sum(alpha) = 0, that is,
    # such that prod(phi) = 1.

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # A fit function must return: theta, sdTheta, thetaOutliers, 
    # phi, sdPhi, phiOutliers.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if (is.null(se)) {
      # For affyPLM v1.10.0 (2006-09-26) or older.
      sdTheta <- rep(1, J);
      sdPhi <- rep(1, I);
    } else {
      # For affyPLM v1.11.6 (2006-11-01) or newer.
      sdTheta <- 2^(se[1:J]);
      sdPhi <- 2^(se[(J+1):length(se)]);
    }
    thetaOutliers <- rep(FALSE, J);
    phiOutliers <- rep(FALSE, I);

    # Return data on the intensity scale
    list(theta=theta, sdTheta=sdTheta, thetaOutliers=thetaOutliers, 
         phi=phi, sdPhi=sdPhi, phiOutliers=phiOutliers);   
  } # rmaModel01()
  attr(rmaModel01, "name") <- "rmaModel01";


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # rmaModel02().
  # Author: Ken Simpson, WEHI, 2006-09-26.
  # Requires: affyPLM() by Ben Bolstad.
  # Why: The above "R_rlm_rma_default_model" call is not available on all
  # platforms (yet).  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rmaModel02 <- function(y, constraint.type=list(default="contr.treatment", chip="contr.treatment", probe="contr.sum")) {
    # Assert right dimensions of 'y'.
    if (length(dim(y)) != 2) {
      str(y);
      stop("Argument 'y' must have two dimensions: ", 
                                                paste(dim(y), collapse="x"));
    }

    # Log-additive model
    y <- log(y, base=2)

    # make factor variables for chip and probe
    nchip <- ncol(y)
    nprobe <- nrow(y)

    chip <- factor(rep(1:nchip, each=nprobe))
    probe <- factor(rep(1:nprobe, nchip))
    X <- model.matrix(~ -1 + chip + probe, contrasts.arg=list(chip=constraint.type$chip, probe=constraint.type$probe))

    # Fit model using affyPLM code
    fit <- .C("rlm_fit_R", as.double(X), as.double(y), rows=as.integer(nchip*nprobe), cols=as.integer(nchip+nprobe-1), beta=double(nchip+nprobe-1), resids=double(nchip*nprobe), weights=double(nchip*nprobe), PACKAGE="affyPLM")

    # Extract probe affinities and chip estimates
    J <- ncol(y);
    I <- nrow(y);
    est <- fit$beta;

    # Chip effects
    beta <- est[1:J];

    # Probe affinities
    alpha <- c(0, est[(J+1):length(est)]);
    if (constraint.type$probe=="contr.sum") {
      alpha[1] <- -sum(alpha[2:length(alpha)]);
    } 
      
    # Estimates on the intensity scale
    theta <- 2^beta;
    phi <- 2^alpha;

    # The RMA model is fitted with constraint sum(alpha) = 0, that is,
    # such that prod(phi) = 1.

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # A fit function must return: theta, sdTheta, thetaOutliers, 
    # phi, sdPhi, phiOutliers.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    sdTheta <- rep(1, J);
    thetaOutliers <- rep(FALSE, J);
    sdPhi <- rep(1, I);
    phiOutliers <- rep(FALSE, I);

    # Return data on the intensity scale
    list(theta=theta, sdTheta=sdTheta, thetaOutliers=thetaOutliers, 
         phi=phi, sdPhi=sdPhi, phiOutliers=phiOutliers);   
  } # rmaModel02()
  attr(rmaModel02, "name") <- "rmaModel02";

  # Test which of the above fit functions are available on the current
  # system, using a dummy set of signals
  rmaModel <- rmaModel02;
  tryCatch({
    rmaModel01(matrix(1:6+0.1, ncol=3));
    rmaModel <- rmaModel01;
  }, error = function(ex) {})

  rmaModel;
}, protected=TRUE)




############################################################################
# HISTORY:
# 2006-11-02
# o Added SE estimates in RmaPlm from Ben's new code. Works with 
#   affyPLM v1.11.6 or newer. /KS
# 2006-09-26
# o Added code to use either of the two RMA fit functions.
# o Incorporated Ken Simpson's fit function for RMA as an alternative.
# 2006-09-11
# o The fit function now returns all required fields.
# 2006-08-25
# o Created from the corresponding Li & Wong model.
############################################################################
