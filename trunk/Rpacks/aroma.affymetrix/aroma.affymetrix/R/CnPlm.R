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
# \arguments{
#   \item{...}{Arguments passed to @see "SnpPlm".}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# \details{
#   Models implementing this copy-number PLM, provides either 
#   allele-specific or total copy-number estimates.  
#   For allele-specific CNs the underlying @see "SnpPlm" model is fitted as
#   is, i.e. for each allele seperately with or without the strands first
#   being merged.
#
#   For total CNs the probe signals for the two alleles are combined 
#   (=summed; not averaged) on the intensity scale before fitting 
#   underlying @see "SnpPlm" model, again with or without the strands 
#   first being merged.
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


setMethodS3("getParameterSet", "CnPlm", function(this, ...) {
  params <- NextMethod("getParameterSet", this, ...);
  params$combineAlleles <- this$combineAlleles;
  params;
}, private=TRUE)


setMethodS3("getCellIndices", "CnPlm", function(this, ...) {
  cells <- NextMethod("getCellIndices", this, ...);

  # If combining alleles, still return all groups as is.
  # The summing is taken care of by the fitUnit() function.
  
  cells;
})



## setMethodS3("getSubname", "CnPlm", function(this, ...) {
##   s <- NextMethod("getSubname", this, ...);
##   if (this$combineAlleles) {
##     s <- sprintf("%sTotal", s);
##   } else {
##     s <- sprintf("%sAllelic", s);
##   }
##   s;
## }, private=TRUE)


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
        yA <- .subset2(.subset2(groups, 1), 1);
        yB <- .subset2(.subset2(groups, 2), 1);
        y <- yA + yB;
        if (length(dim(y)) == 3) {
          y <- y[1,,] + y[2,,];
        }
        list(fitfcn(y));
      } else if (ngroups == 4) {
        yA1 <- .subset2(.subset2(groups, 1), 1);
        yB1 <- .subset2(.subset2(groups, 2), 1);
        yA2 <- .subset2(.subset2(groups, 3), 1);
        yB2 <- .subset2(.subset2(groups, 4), 1);
        y1 <- yA1 + yB1;
        y2 <- yA2 + yB2;
        if (length(dim(y1)) == 3) {
          y1 <- y1[1,,] + y1[2,,];
        }
        if (length(dim(y2)) == 3) {
          y2 <- y2[1,,] + y2[2,,];
        }
        list(
          fitfcn(y1), 
          fitfcn(y2)
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
}, private=TRUE)


setMethodS3("getChipEffectSetClass", "CnPlm", function(this, ...) {
  CnChipEffectSet;
}, private=TRUE)


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
# 2006-12-10
# o Added support to fit PLM to MM or (PM+MM).
# 2006-09-11
# o Recreated.
############################################################################
