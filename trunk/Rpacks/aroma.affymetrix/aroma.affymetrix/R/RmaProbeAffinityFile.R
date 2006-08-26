###########################################################################/**
# @RdocClass RmaProbeAffinityFile
#
# @title "The RmaProbeAffinityFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of probe affinities in the 
#  RMA model, see @see "AffymetrixRmaModel".
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
#   This class subclasses the @see "ProbeAffinityFile" class so that the 
#   above estimates can be stored in a CEL file.  Please see that class for 
#   more details on the general idea behind encoding/decoding data in
#   CEL files.
#
#   The mapping/encoding is as follows:
#   \itemize{
#     \item{intensities}{Hold \code{phi} (probe affinities).}
#     \item{stdvs}{Hold \code{stdvs} (SD of the probe affinities).}
#     \item{pixels}{Currently just zeros.}
#   }
#   For more details about the file format, see also the source code of 
#   this class.
# }
# 
# @author
#*/###########################################################################
setConstructorS3("RmaProbeAffinityFile", function(...) {
  extend(ProbeAffinityFile(...), "RmaProbeAffinityFile")
})


setMethodS3("encodeUnitGroup", "RmaProbeAffinityFile", function(static, groupData, ...) {
  phi <- .subset2(groupData, "phi");
  ncells <- length(phi);
  pixels <- rep(0, ncells);
  stdvs <- rep(1, ncells);
  list(intensities=phi, stdvs=stdvs, pixels=pixels);
}, static=TRUE, protected=TRUE)




setMethodS3("decodeUnitGroup", "RmaProbeAffinityFile", function(static, groupData, ...) {
  attachLocally(groupData);

  pixels <- groupData$pixels;

  list(phi=groupData$intensities, stdvs=groupData$stdvs, pixels=pixels);
}, static=TRUE, protected=TRUE)


############################################################################
# HISTORY:
# 2006-08-25
# o Created.
############################################################################
