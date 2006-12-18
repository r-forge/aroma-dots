###########################################################################/**
# @RdocClass MbeiPlm
#
# @title "The MbeiPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Li \& Wong (2001) model.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeLevelModel".}
#   \item{tags}{A @character @vector of tags.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#   For a single unit group, the multiplicative model of dChip is:
#
#    \deqn{y_{ij} = \theta_i \phi_j + \varepsilon_{ij}}
#
#   where \eqn{\theta_i} are the chip effects for arrays \eqn{i=1,...,I}, 
#   and \eqn{\phi_j} are the probe affinities for probes \eqn{j=1,...,J}.
#   The \eqn{\varepsilon_{ij}} are zero-mean noise with equal variance.
#
#   In addition, we modify th constraint such that it is guaranteed that
#   \eqn{\prod_j \phi_j = 1}.
# }
#
# @author
#
# \references{
#   Li, C. and Wong, W.H. (2001), Genome Biology 2, 1-11.\cr
#   Li, C. and Wong, W.H. (2001), Proc. Natl. Acad. Sci USA 98, 31-36.\cr
# }
#*/###########################################################################
setConstructorS3("MbeiPlm", function(..., tags="*") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(affy) || throw("Package 'affy' not loaded.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    tags[tags == "*"] <- "MBEI";
  }


  extend(ProbeLevelModel(..., tags=tags), "MbeiPlm")
})



setMethodS3("getProbeAffinities", "MbeiPlm", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the probe affinities (and create files etc)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  paf <- NextMethod("getProbeAffinities", this, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the encode and decode functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setEncodeFunction(paf, function(groupData, ...) {
    # Rename some fields so that we support the structure of this class,
    # but also output from affy::fit.li.wong().
    names <- names(groupData);
    # Is it an affy:fit.li.wong() structure?
    if ("sdPhi" %in% names) {
      names <- sub("iter", "nbrOfIterations", names);
      names <- sub("convergence1", "converged", names);
      names <- sub("convergence2", "convergedOutliers", names);
      names(groupData) <- names;
    }
  
    # Encode outliers as the sign of 'pixels'; -1 = TRUE, +1 = FALSE
    pixels <- sign(0.5 - as.integer(groupData$phiOutliers));
      
    # Encode the number of iterations as the absolute value of the 1st pixel.
    pixels[1] <- pixels[1]*groupData$nbrOfIterations;
      
    # Encode convergence1/2 as bits in the 2nd pixel.
    # Note: For this to work, there must be at least two 'pixels' in the
    # unit group.  However, this is not true for some strands on the 500K
    # chip, which might have a single probe quartet, e.g. SNP_A-1781633 on
    # the Nsp chip.  Thus, this is a problem if the MBEI model is fitted
    # strand by strand.  On the other hand, the PLM based on chip effects
    # and probe affinities breaks down if there is only one probe, so 
    # talking about probe affinities does not make sense in the first 
    # place.  Thus, the probe affinity for this single probe will be set
    # to exactly one (on the intensity scale) by the construct that the 
    # probe affinities should multiply to one.  There is no standard 
    # deviation estimate or convergence results for such single-probe 
    # models. /HB 2006-12-18
    npixels <- length(pixels);
    if (npixels > 1) {
      pixels[2] <- pixels[2] * 
              (1 + 2*groupData$converged + 4*groupData$convergedOutliers);
    }
  
    list(intensities=groupData$phi, stdvs=groupData$sdPhi, pixels=pixels);
  })

  setDecodeFunction(paf,  function(groupData, ...) {
    pixels <- groupData$pixels;
  
    # Outliers are encoded by the sign of 'pixels'.
    outliers <- as.logical(1-sign(pixels));
  
    # Number of iterations is encoded as the absolute value of the 1st pixel.
    nbrOfIterations <- as.integer(abs(pixels[1])+0.5);

    npixels <- length(pixels);
    if (npixels > 1) {
      # convergence & convergenceOutliers are encoded as bits in the 2nd pixel.
      t <- pixels[2] %/% 2;
      converged <- as.logical(t %% 2 == 1);  t <- t %/% 2;
      convergedOutliers <- as.logical(t %% 2 == 1);
    } else {
      # See comments on the encode function above. /HB 2006-12-18
      converged <- TRUE;
      convergedOutliers <- TRUE;
    }
  
    list(
      phi=groupData$intensities, 
      sdPhi=groupData$stdvs, 
      phiOutliers=outliers, 
      nbrOfIterations=nbrOfIterations, 
      converged=converged, 
      convergedOutliers=convergedOutliers
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
#   @see "affy::fit.li.wong".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitFunction", "MbeiPlm", function(this, ...) {
  standardize <- this$standardize;

  liWong <- function(y, ...) {
    # Enough of probes?
    if (nrow(y) > 1) {
      fit <- fit.li.wong(t(y));

      # A fit function must return: theta, sdTheta, thetaOutliers, 
      # phi, sdPhi, phiOutliers.
      names <- names(fit);
      idxs <- match(c("sigma.theta", "theta.outliers", "sigma.phi", 
                                                     "phi.outliers"), names);
      names[idxs] <- c("sdTheta", "thetaOutliers", "sdPhi", "phiOutliers");
      names(fit) <- names;
  
      # Rescale such that prod(phi) = 1?
      if (standardize) {
        phi <- fit$phi;
        theta <- fit$theta;
        I <- length(phi);
        c <- prod(phi)^(1/I);
        phi <- phi/c;
        theta <- theta*c;
        fit$phi <- phi;
        fit$theta <- theta;
      }
    } else {
      # For the case where there is only a single probe in the unit group
      # let the chip effect be the probe signal and the probe affinity one.
      # This is indeed what fit.li.wong() returns, but we don't want their
      # warning about it:
      fit <- list(theta=as.vector(y), sdTheta=NA, thetaOutliers=NA, 
                  phi=1, sdPhi=NA, phiOutliers=NA, 
                  sigma.eps=NA, single.outliers=NA,
                  convergence1=NA, convergence2=NA, 
                  iter=1, delta=NA);
    }

    fit;
  }

  liWong;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2006-12-18
# o Now the fit function of the MBEI model treats single-probe unit groups
#   specially; the affy::fit.li.wong() handled it already before, but 
#   generated a warning for each call.
# o BUG FIX: The encoding/decoding of probe affinities assumed at least two
#   probes per unit group, but this is not true for all chip types, e.g.
#   the 500K SNP chips.
# 2006-09-10
# o Updated getFitFunction() to return required fields.
# 2006-08-24
# o Added Rdoc comments.
# 2006-08-23
# o Added getProbeAffinities() and the corrsponding cached fields.
# o Now fit() does not re-read data just updated.
# 2006-08-19
# o After all the bug fixes in updateCel() I think this function finally
#   works correctly.
# 2006-08-17
# o Created.
############################################################################
