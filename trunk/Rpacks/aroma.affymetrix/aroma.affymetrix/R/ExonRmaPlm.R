###########################################################################/**
# @RdocClass ExonRmaPlm
#
# @title "The ExonRmaPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the log-additive model part of the Robust Multichip
#  Analysis (RMA) method described in Irizarry et al (2003), as implemented
#  for exon arrays.  The model may be fitted with exons merged into
#  transcripts (all probes fitted together) or on an individual exon basis
#  (probes within an exon treated as a group, but exons fitted separately).
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RmaPlm".}
#   \item{tags}{A @character @vector of tags.}
#   \item{mergeGroups}{A @logical flag specifying whether to merge exons
#      into transcripts.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#    @see "RmaPlm".
# }
#
# \author{Ken Simpson (ksimpson[at]wehi.edu.au).}
#
# \references{
#  Irizarry et al. \emph{Summaries of Affymetrix GeneChip probe level data}. 
#  NAR, 2003, 31, e15.\cr
# }
#*/###########################################################################
setConstructorS3("ExonRmaPlm", function(..., tags="*", mergeGroups=TRUE) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    asteriskTag <- "RMA";
    # Update default tags
    tags[tags == "*"] <- asteriskTag;

    if (mergeGroups) {
      tags <- c(tags, "merged");
    }
    
    # Split by commas
    tags <- paste(tags, collapse=",");
    tags <- unlist(strsplit(tags, split=","));
  }

  extend(RmaPlm(..., tags=tags), "ExonRmaPlm",
    mergeGroups=mergeGroups
  )
})

# utility function - keep here for now

cdfMergeGroups <- function(groups, ...) {
  nbrOfGroups <- length(groups);
  
  nbrOfFields <- length(.subset2(groups,1));
  newGroup <- vector("list", nbrOfFields);
  for (ff in seq(length=nbrOfFields)) {
    newGroup[[ff]] <- unlist(base::lapply(groups, .subset2, ff), use.names=FALSE);
  }
  names(newGroup) <- names(.subset2(groups,1));
  return(list(newGroup));
}

setMethodS3("getCellIndices", "ExonRmaPlm", function(this, ...) {

  cells <- NextMethod("getCellIndices", this, ...);

  # Merge groups?
  if (this$mergeGroups) {
    cells <- applyCdfGroups(cells, cdfMergeGroups);
  }

  cells;
})


setMethodS3("getChipEffectSetClass", "ExonRmaPlm", function(this, ...) {
  ExonChipEffectSet;
}, private=TRUE)


setMethodS3("getChipEffectSet", "ExonRmaPlm", function(this, ...) {
  ces <- NextMethod("getChipEffectSet", this, ...);
  setMergeGroups(ces, this$mergeGroups);
  ces;
})

setMethodS3("getChipEffects", "ExonRmaPlm", function(this, ...) {
  getChipEffectSet(this, ...);
})



setMethodS3("getProbeAffinityFile", "ExonRmaPlm", function(this, ..., .class=ExonProbeAffinityFile) {
  paf <- NextMethod("getProbeAffinityFile", this, ..., .class=.class);
  setMergeGroups(paf, this$mergeGroups);
  paf;
})


setMethodS3("getProbeAffinities", "ExonRmaPlm", function(this, ...) {
  getProbeAffinityFile(this, ...);
})


setMethodS3("setMergeGroups", "ExonRmaPlm", function(this, ...) {
  ces <- getChipEffects(this);
  setMergeGroups(ces, ...);
  paf <- getProbeAffinities(this);
  setMergeGroups(paf, ...);
})


###########################################################################/**
# @RdocMethod getFitFunction
#
# @title "Gets the low-level function that fits the Exon PLM"
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
# \author{
#   Henrik Bengtsson and Ken Simpson (WEHI) utilizing Ben Bolstad's 
#   \pkg{affyPLM} package.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitFunction", "ExonRmaPlm", function(this, ...) {
  # This should not be need, but for some reason is the package not loaded
  # although it is listed in DESCRIPTION. /HB 2007-02-09
  require("affyPLM") || throw("Package not loaded: affyPLM");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # rmaModelAffyPlm()
  # Author: Henrik Bengtsson, UC Berkeley. 
  # Requires: affyPLM() by Ben Bolstad.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  exonRmaModel <- function(y, psiCode=0, psiK=1.345, useMedianPolish=FALSE, medianPolishThreshold=100){
    # Assert right dimensions of 'y'.

    # If input data are dimensionless, return NAs. /KS 2006-01-30
    if (is.null(dim(y))) {
      nbrOfArrays <- nbrOfArrays(getDataSet(this));
      return(list(theta=rep(NA, nbrOfArrays),
                  sdTheta=rep(NA, nbrOfArrays),
                  thetaOutliers=rep(NA, nbrOfArrays), 
                  phi=c(), 
                  sdPhi=c(), 
                  phiOutliers=c()));
      
    }

    if (length(dim(y)) != 2) {
      str(y);
      stop("Argument 'y' must have two dimensions: ", 
                                                paste(dim(y), collapse="x"));
    }

    # Log-additive model
    y <- log(y, base=2);

    # do we want to use median polish for large matrices?

    if (useMedianPolish) {
      nbrOfVariables <- dim(y)[1]+dim(y)[2]-1;
      if (nbrOfVariables > medianPolishThreshold) {
        mp <- medpolish(y, trace.iter=FALSE);
        fit <- list(Estimates=c(mp$overall+mp$col, mp$row), StdErrors=rep(0, length(c(mp$row,mp$col))));
      } else {    
      # Fit model using affyPLM code
        fit <- .Call("R_rlm_rma_default_model", y, psiCode, psiK, PACKAGE="affyPLM");
      }
    } else {
      fit <- .Call("R_rlm_rma_default_model", y, psiCode, psiK, PACKAGE="affyPLM");      
    }
        
    # Extract probe affinities and chip estimates
    J <- ncol(y);  # Number of arrays
    I <- nrow(y);  # Number of probes
    est <- fit$Estimates;
    se <- fit$StdErrors;

    # Chip effects
    beta <- est[1:J];

    # Probe affinities.  If only one probe, must have affinity=1 since
    # sum constraint => affinities sum to zero (on log scale)
    if (I==1) {
      alpha <- 0;
    } else {
      alpha <- est[(J+1):length(est)];
      alpha[length(alpha)] <- -sum(alpha[1:(length(alpha)-1)]);
    }
      
    # Estimates on the intensity scale
    theta <- 2^beta;
    phi <- 2^alpha;

    # The RMA model is fitted with constraint sum(alpha) = 0, that is,
    # such that prod(phi) = 1.

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # A fit function must return: theta, sdTheta, thetaOutliers, 
    # phi, sdPhi, phiOutliers.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.null(se)) {
      # For affyPLM v1.10.0 (2006-09-26) or older.
      sdTheta <- rep(1, J);
      sdPhi <- rep(1, I);
    } else {
      # For affyPLM v1.11.6 (2006-11-01) or newer.
      sdTheta <- 2^(se[1:J]);
      sdPhi <- 2^(se[(J+1):length(se)]);
    }
    thetaOutliers <- rep(FALSE, J);
    phiOutliers <- rep(FALSE, I);

    # Return data on the intensity scale
    list(theta=theta, sdTheta=sdTheta, thetaOutliers=thetaOutliers, 
         phi=phi, sdPhi=sdPhi, phiOutliers=phiOutliers);   
  } # rmaModelAffyPlm()
  attr(exonRmaModel, "name") <- "exonRmaModel";

  exonRmaModel;
}, private=TRUE)




##############################################################################
# HISTORY:
# 2007-07-13 /HB
# o Removed findUnitsTodo() of ExonRmaPlm; it just replicated the superclass.
# 2007-04-24 /HB+EP
# o getProbeAffinityFile() of ExonRmaPlm did not return the correct subclass.
# 2007-04-13 /HB+EP
# o BUG FIX: getChipEffectSet() and getProbeAffinityFile() did not set the
#   'mergeStrands' parameter.  Thanks Elizabeth Purdom for the fix.
# 2006-??-??
# o Created by Ken Simpson, WEHI.
##############################################################################
