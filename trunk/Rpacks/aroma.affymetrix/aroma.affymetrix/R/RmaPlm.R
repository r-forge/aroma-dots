###########################################################################/**
# @RdocClass RmaPlm
#
# @title "The RmaPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the log-additive model used in RMA.
#  It can be used to fit the model on a @see "AffymetrixCelSet".
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeLevelModel".}
#   \item{name}{The name of the model, which is also used in the pathname.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#   Consider a specific unit group.  The log-additive model in RMA is:
#
#    \deqn{log_2(y_{ij})} = \beta_i + \alpha_j + \eps_{ij}
#
#   where \eqn{\beta_i} are the chip effects for arrays \eqn{i=1,...,I}, 
#   and \eqn{\alpha_j} are the probe affinities for probes \eqn{j=1,...,J}.
#   The \eqn{\eps_{ij}} are zero-mean noise with equal variance.
#
#   To minimize the risk for mistakes, all probe-affinity models return
#   parameter estimates on the intensity scale.  That is, this class
#   returns \eqn{\theta_i = 2^\beta_i} and \eqn{\phi_i = 2^\alpha_i},
#   cf. the multiplicative model of Li & Wong.
#
#   Use @seemethod "getProbeAffinities" to get the probe-affinity estimates.
# }
#
# @author
#
# \section{Model estimates}{
#   The estimated probe affinities are represented by the
#   @see "RmaProbeAffinityFile" class.  
# }
#
# \references{
# }
#*/###########################################################################
setConstructorS3("RmaPlm", function(..., name="modelRmaPlm") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(affyPLM) || throw("Package 'affyPLM' not loaded.");

  extend(ProbeLevelModel(..., name=name), "RmaPlm")
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
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitFunction", "RmaPlm", function(this, ...) {
  rmaModel <- function(y, psiCode=0, psiK=1.345){
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
  }

  rmaModel;
}, protected=TRUE)




############################################################################
# HISTORY:
# 2006-09-11
# o The fit function now returns all required fields.
# 2006-08-25
# o Created from the corresponding Li & Wong model.
############################################################################
