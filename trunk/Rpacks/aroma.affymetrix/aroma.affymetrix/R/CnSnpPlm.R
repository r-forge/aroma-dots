###########################################################################/**
# @RdocClass CnSnpPlm
#
# @title "The CnSnpPlm class"
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
#   For total CNs the probe signals for the two alleles are averaged on the
#   intensity scale before fitting underlying @see "SnpPlm" model, again
#   with or without the strands first being merged.
# }
#
# \section{Requirments}{
#   Classes inheriting from this @see "Interface" must provide the following
#   fields, in addition to the ones according to @see "SnpPlm":
#   \itemize{
#    \item{averageAB}{A @logical indicating if total or allele-specific
#      copy numbers should be estimated according to the above averaging.}
#   }
# }
#
# @author
#*/###########################################################################
setConstructorS3("CnSnpPlm", function(...) {
  extend(SnpPlm(...), "CnSnpPlm");
})


setMethodS3("getFitUnitFunction", "CnSnpPlm", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Select fit function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Total copy number of not?
  if (this$averageAB) {
    # Get the fit function for a single set of intensities
    fitfcn <- getFitFunction(this, ...);
  
    fitUnit <- function(groups, ...) {
      ngroups <- length(groups);
      if (ngroups == 2) {
        yA <- .subset2(groups, 1);
        yB <- .subset2(groups, 2);
        list(fitfcn((yA + yB)/2));
      } else if (ngroups == 4) {
        yA1 <- .subset2(groups, 1);
        yB1 <- .subset2(groups, 2);
        yA2 <- .subset2(groups, 3);
        yB2 <- .subset2(groups, 4);
        list(
          fitfcn((yA1 + yB1)/2), 
          fitfcn((yA2 + yB2)/2)
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

setMethodS3("getChipEffects", "CnSnpPlm", function(this, ...) {
  ces <- NextMethod("getChipEffects", this, ...);
  setAverageAB(ces, this$averageAB);
  ces;
})

setMethodS3("getProbeAffinities", "CnSnpPlm", function(this, ..., .class=CnSnpProbeAffinityFile) {
  paf <- NextMethod("getProbeAffinities", this, ..., .class=.class);
  setAverageAB(paf, this$averageAB);
  paf;
})

setMethodS3("setAverageAB", "CnSnpPlm", function(this, ...) {
  ces <- getChipEffects(this);
  setAverageAB(ces, ...);
  paf <- getProbeAffinities(this);
  setAverageAB(paf, ...);
})


############################################################################
# HISTORY:
# 2006-09-11
# o Recreated.
############################################################################
