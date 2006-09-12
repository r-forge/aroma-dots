###########################################################################/**
# @RdocClass CnPlm
#
# @title "The CnPlm class"
#
# \description{
#  @classhierarchy
#
#  This suport class represents a @see "SnpPlm" specially designed for 
#  copy-number analysis.
# }
# 
# @synopsis
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# \details{
#   Models implementing this copy-number PLM, provides either allele-specific
#   or total copy-number estimates.  
#   For allele-specific CNs the underlying @see "SnpPlm" model is fitted as
#   is, i.e. for each allele seperately with or without the strands first
#   being merged.
#
#   For total CNs the probe signals for the two alleles are combined 
#   (=summed) on the intensity scale before fitting underlying 
#   @see "SnpPlm" model, again with or without the strands first being
#   merged.
# }
#
# \section{Requirments}{
#   Classes inheriting from this @see "Interface" must provide the following
#   fields, in addition to the ones according to @see "SnpPlm":
#   \itemize{
#    \item{combineAlleles}{A @logical indicating if total or allele-specific
#      copy numbers should be estimated according to the above averaging.}
#   }
# }
#
# @author
#*/###########################################################################
setConstructorS3("CnPlm", function(...) {
  extend(SnpPlm(...), "CnPlm");
})


setMethodS3("getFitUnitFunction", "CnPlm", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Select fit function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Total copy number of not?
  if (this$combineAlleles) {
    # Get the fit function for a single set of intensities
    fitfcn <- getFitFunction(this, ...);
  
    fitUnit <- function(groups, ...) {
      ngroups <- length(groups);
      if (ngroups == 2) {
        yA <- .subset2(groups, 1);
        yB <- .subset2(groups, 2);
        list(fitfcn(yA + yB));
      } else if (ngroups == 4) {
        yA1 <- .subset2(groups, 1);
        yB1 <- .subset2(groups, 2);
        yA2 <- .subset2(groups, 3);
        yB2 <- .subset2(groups, 4);
        list(
          fitfcn(yA1 + yB1), 
          fitfcn(yA2 + yB2)
        );
      } else {
        # For all other cases, fit each group individually
        lapply(groups, FUN=function(group) {
          y <- .subset2(group, 1);
          fitfcn(y);
        })
      }
    }
  } else {
    fitUnit <- NextMethod("getFitUnitFunction", this, ...);
  }

  fitUnit;
})

setMethodS3("getChipEffectSetClass", "CnPlm", function(this, ...) {
  CnChipEffectSet;
})

setMethodS3("getChipEffects", "CnPlm", function(this, ...) {
  ces <- NextMethod("getChipEffects", this, ...);
  setCombineAlleles(ces, this$combineAlleles);
  ces;
})

setMethodS3("getProbeAffinities", "CnPlm", function(this, ..., .class=CnProbeAffinityFile) {
  paf <- NextMethod("getProbeAffinities", this, ..., .class=.class);
  setCombineAlleles(paf, this$combineAlleles);
  paf;
})

setMethodS3("setCombineAlleles", "CnPlm", function(this, ...) {
  ces <- getChipEffects(this);
  setCombineAlleles(ces, ...);
  paf <- getProbeAffinities(this);
  setCombineAlleles(paf, ...);
})


############################################################################
# HISTORY:
# 2006-09-11
# o Recreated.
############################################################################
