###########################################################################/**
# @RdocClass AffinePlm
#
# @title "The AffinePlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Bengtsson \& Hössjer (2006) model,
#  see @see "AffinePlm".
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeLevelModel".}
#   \item{name}{The name of the model, which is also used in the pathname.}
#   \item{background}{If @TRUE, background is estimate for each unit group,
#     otherwise not. That is, if @FALSE, a \emph{linear} model without 
#     offset is fitted, resulting in very similar results as obtained by
#     the @see "MbeiPlm".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \references{
#   Bengtsson \& Hössjer (2006). \cr
# }
#*/###########################################################################
setConstructorS3("AffinePlm", function(..., name="modelAffinePlm", background=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(aroma.light) || throw("Package 'aroma.light' not loaded.");

  extend(ProbeLevelModel(..., name=name), "AffinePlm",
    background = background
  )
})


setMethodS3("getProbeAffinities", "AffinePlm", function(this, ...) {
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
#   @see "aroma.light::calibrateMultiscan".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitFunction", "AffinePlm", function(this, ...) {
  standardize <- this$standardize;
  center <- this$background;

  affineFit <- function(y, ...) {
    # NOTE: If center=FALSE => constraint a=0 /HB 2006-09-11
    f <- calibrateMultiscan(t(y), center=center, project=TRUE);
    theta <- as.vector(f);
    phi <- as.vector(attr(f, "modelFit")$b);

    J <- length(theta);
    I <- length(phi);

    # Rescale such that prod(phi) = 1?
    if (standardize) {
      c <- prod(phi)^(1/I);
      phi <- phi/c;
      theta <- theta*c;
    }

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
  } # affineFit()

  affineFit;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2006-09-11
# o Added argument 'background' to fit background or not.
# 2006-08-28
# o Created from the Li & Wong model.
############################################################################
