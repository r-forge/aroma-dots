###########################################################################/**
# @RdocClass AffinePlm
#
# @title "The AffinePlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Bengtsson \& Hössjer (2006) model.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeLevelModel".}
#   \item{background}{If @TRUE, background is estimate for each unit group,
#     otherwise not. That is, if @FALSE, a \emph{linear} model without 
#     offset is fitted, resulting in very similar results as obtained by
#     the @see "MbeiPlm".}
#   \item{tags}{A @character @vector of tags.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#   For a single unit group, the affine model is:
#
#    \deqn{y_{ij} = a + \theta_i \phi_j + \varepsilon_{ij}}
#
#   where \eqn{a} is an offset common to all probe signals, 
#   \eqn{\theta_i} are the chip effects for arrays \eqn{i=1,...,I}, 
#   and \eqn{\phi_j} are the probe affinities for probes \eqn{j=1,...,J}.
#   The \eqn{\varepsilon_{ij}} are zero-mean noise with equal variance.
#   The model is constrained such that \eqn{\prod_j \phi_j = 1}.
#
#   Note that with the additional constraint \eqn{a=0} (see arguments above),
#   the above model is very similar to @see "MbeiPlm".  The differences in
#   parameter estimates is due to difference is assumptions about the
#   error structure, which in turn affects how the model is estimated.
# }
#
# @author
#
# \references{
#   Bengtsson \& Hössjer (2006). \cr
# }
#*/###########################################################################
setConstructorS3("AffinePlm", function(..., background=TRUE, tags="*") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(aroma.light) || throw("Package 'aroma.light' not loaded.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    idx <- which(tags == "*");
    if (length(idx) > 0) {
      tags[idx] <- "APLM";
      if (!background)
        tags <- R.utils::insert.default(tags, idx+1, "linear");
    }
  }


  extend(ProbeLevelModel(..., tags=tags), "AffinePlm",
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
