###########################################################################/**
# @RdocClass AffymetrixRmaModel
#
# @title "The AffymetrixRmaModel class"
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
#   \item{dataSet}{An @see "AffymetrixCelSet" object.}
#   \item{path}{The @character string specifying the path to the directory
#      to contain the parameter-estimate files.}
#   \item{...}{Arguments passed to @see "ProbeLevelModel".}
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
setConstructorS3("AffymetrixRmaModel", function(dataSet=NULL, path=filePath("modelRMA", getChipType(getCdf(dataSet))), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(affyPLM) || throw("Package 'affyPLM' not loaded.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  if (is.null(dataSet)) {
    # A work-around for the fact that getCdf(NULL) is not work.
    path=NULL;
  }

  extend(ProbeLevelModel(dataSet=dataSet, path=path, ...), "AffymetrixRmaModel")
})


setMethodS3("getProbeAffinityClass", "AffymetrixRmaModel", function(static, ...) {
  RmaProbeAffinityFile;
}, static=TRUE, protected=TRUE)



setMethodS3("getFitFunction", "AffymetrixRmaModel", function(static, ...) {
  rmaModel <- function(y, psiCode=0, psiK=1.345){
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

#    list(beta=beta, alpha=alpha);

    # Return data on the intensity scale
    list(theta=theta, phi=phi);   
  }

  rmaModel;
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2006-08-25
# o Created from the corresponding Li & Wong model.
############################################################################
