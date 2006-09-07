###########################################################################/**
# @RdocClass RmaSnpPlm
#
# @title "The RmaSnpPlm class"
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
setConstructorS3("RmaSnpPlm", function(..., name="rmaSnpPlm") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(affyPLM) || throw("Package 'affyPLM' not loaded.");

  extend(RmaPlm(..., name=name), "RmaSnpPlm")
})



setMethodS3("getFitFunction", "RmaSnpPlm", function(static, ...) {
  rmaModel <- function(y, psiCode=0, psiK=1.345){
    # Assert right dimensions of 'y'.
    if (length(dim(y)) != 2) {
      str(y);
      stop("Argument 'y' must have two dimensions: ", paste(dim(y), collapse="x"));
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
