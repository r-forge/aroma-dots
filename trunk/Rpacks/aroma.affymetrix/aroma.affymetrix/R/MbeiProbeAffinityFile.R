###########################################################################/**
# @RdocClass MbeiProbeAffinityFile
#
# @title "The MbeiProbeAffinityFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of probe affinities in the 
#  Li \& Wong (2001) model.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ParameterCelFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Stored parameters}{
#   For any given unit group consisting of L probes (L >= 2), the
#   following information is stored:
#   \itemize{
#     \item{phi}{(L @doubles [floats] in (0,Inf)) The probe affinities, 
#       one for each probe.}
#     \item{stdvs}{(L @doubles [floats] in (0,Inf)) The standard deviation
#       of the probe-affinity estimates, one for each probe.}
#     \item{outliers}{(L logicals) Specifies if a probe was considered to 
#       be an outlier or not, one for each probe.}
#     \item{nbrOfIterations}{(an @integer in [0,Inf]) Number of iterations
#       before the algorithm converged.  
#       A \code{0} indicates that this unit group has not been fitted.}
#     \item{converged}{(a @logical) If @FALSE, the algorithm did not 
#       converge when fitting the phis (and thetas).}
#     \item{convergedOutliers}{(a @logical) If @FALSE, the algorithm 
#       did not converge in deciding what are outliers.}
#   }
#
#   For more details about these parameters, see @see "affy::fit.li.wong".
# }
#
# \section{File format}{
#   This class subclasses @see "ParameterCelFile" so that the above 
#   estimates can be stored in a CEL file.  Please see that class for 
#   more details on the general idea behind encoding/decoding data in
#   CEL files.
#
#   The mapping/encoding is as follows:
#   \itemize{
#     \item{intensities}{Hold \code{phi} (probe affinities).}
#     \item{stdvs}{Hold \code{stdvs} (SD of the probe affinities).}
#     \item{pixels}{Hold \code{outliers}, \code{nbrOfIterations}, 
#       \code{converged}, and \code{convergedOutliers}.
#       This is done by letting the sign code for \code{outliers} such that
#       a negative sign corresponds to @TRUE and a positive to @FALSE.
#       Then the absolute value of \code{pixels[1]} equals 
#       \code{nbrOfIterations}.
#       The \code{converged} is encoded as the second bit (module 2) of
#       \code{pixel[2]}, and \code{convergedOutliers} as the third 
#       bit (modulo 4).  Note that we know that there are at least two
#       probes, because otherwise we cannot fit the model.}
#   }
#   For more details about the file format, see also the source code of 
#   this class.
# }
# 
# @author
# 
# \seealso{
#   An object of this class is typically obtain through the
#   \code{getProbeAffinities()} method for the 
#   @see "ChipPlm" class.
# }
#
#*/###########################################################################
setConstructorS3("MbeiProbeAffinityFile", function(...) {
  extend(ProbeAffinityFile(...), "MbeiProbeAffinityFile");
})


setMethodS3("encodeUnitGroup", "MbeiProbeAffinityFile", function(static, groupData, ...) {
  # Rename some fields (so that we support the structure of this class,
  # but also output from affy::fit.li.wong().
  names <- names(groupData);
  # Is it an affy:fit.li.wong() structure?
  if ("sigma.phi" %in% names) {
    names <- sub("sigma.phi", "stdvs", names);
    names <- sub("phi.outliers", "outliers", names);
    names <- sub("iter", "nbrOfIterations", names);
    names <- sub("convergence1", "converged", names);
    names <- sub("convergence2", "convergedOutliers", names);
    names(groupData) <- names;
  }

  # Encode outliers as the sign of 'pixels'; -1 = TRUE, +1 = FALSE
  pixels <- sign(0.5 - as.integer(groupData$outliers));
    
  # Note: There are at least two 'pixels', otherwise we can't fit the model.
    
  # Ecode the number of iterations as the absolute value of the 1st pixel.
  pixels[1] <- pixels[1]*groupData$nbrOfIterations;
    
  # Ecode convergence1/2 as bits in the 2nd pixel.
  pixels[2] <- pixels[2] * 
            (1 + 2*groupData$converged + 4*groupData$convergedOutliers);

  list(intensities=groupData$phi, stdvs=groupData$stdvs, pixels=pixels);
}, static=TRUE, protected=TRUE)




setMethodS3("decodeUnitGroup", "MbeiProbeAffinityFile", function(static, groupData, ...) {
  attachLocally(groupData);

str(groupData);

  pixels <- groupData$pixels;

  # Outliers are encoded by the sign of 'pixels'.
  outliers <- as.logical(1-sign(pixels));

  # Number of iterations is encoded as the absolute value of the 1st pixel.
  nbrOfIterations <- as.integer(abs(pixels[1])+0.5);

  # convergence & convergenceOutliers are encoded as bits in the 2nd pixel.
  t <- pixels[2] %/% 2;
  converged <- as.logical(t %% 2 == 1);  t <- t %/% 2;
  convergedOutliers <- as.logical(t %% 2 == 1);

  list(phi=groupData$intensities, stdvs=groupData$stdvs, outliers=outliers, nbrOfIterations=nbrOfIterations, converged=converged, convergedOutliers=convergedOutliers);
}, static=TRUE, protected=TRUE)


############################################################################
# HISTORY:
# 2006-08-24
# o Added details Rdoc comments.
# 2006-08-23
# o Added updateUnits().
# 2006-08-21
# o Created.
############################################################################
