###########################################################################/**
# @RdocClass RmaCnSnpPlm
#
# @title "The RmaCnSnpPlm class"
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
# \section{Model}{
#   Consider a specific unit group.  The log-additive model in RMA is:
#
#    \deqn{log_2(y_{ij})} = \beta_i + \alpha_j + \eps_{ij}
#
#   where \eqn{\beta_i} are the chip effects for arrays \eqn{i=1,...,I}, 
#   and \eqn{\alpha_j} are the probe affinities for probes \eqn{j=1,...,J}.
#   The \eqn{\eps_{ij}} are zero-mean noise with equal variance.
#
#   To minimize the risk for mistakes, all probe-affinity models return
#   parameter estimates on the intensity scale.  That is, this class
#   returns \eqn{\theta_i = 2^\beta_i} and \eqn{\phi_i = 2^\alpha_i},
#   cf. the multiplicative model of Li & Wong.
#
#   Use @seemethod "getProbeAffinities" to get the probe-affinity estimates.
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
setConstructorS3("RmaCnSnpPlm", function(..., name="modelTotalCnRma") {
  this <- extend(RmaPlm(..., name=name), "RmaCnSnpPlm")
  this <- setup(this);
  this;
})


setMethodS3("setup", "RmaCnSnpPlm", function(this, ...) {
  dataSet <- getDataSet(this);
  if (is.null(dataSet))
    return(invisible(this));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdfAB <- function(units) {
    if (is.list(units[[1]])) {
      applyCdfGroups(units, function(groups) {
        groups <- cdfMergeStrands(groups);
        if (length(groups) >= 2)
          groups <- cdfMergeAlleles(groups)
        groups;
      });
    } else {
      lapply(units, FUN=function(n) {
        paste(unique(n), collapse="");
      });
    }
  }
  
  if (!is.null(dataSet)) {
    # Enforce that the right CDF structure is returned (for now at least).
    cdf <- getCdf(dataSet);
    setRestructor(cdf, cdfAB);
    setCdf(dataSet, cdf);
  }

  invisible(this);
}, protected=TRUE);



setMethodS3("fit", "RmaCnSnpPlm", function(this, ..., transform=NULL, postCdfTransform=NULL) {
  # The CDF object should be such that it returns an array where the first
  # dimension has alleles "A" and "B".  Thus, when retrieving the data,
  # data[1,,] is signals for allele A, and data[2,,] for allele B.
  if (is.null(transform)) {
    transform <- function(x) {
      if (length(dim(x)) > 2)
        x <- colMeans(x);
      x;
    }
  }

  # When storing the data we only store one estimate per allele (A,B) pair.
  if (is.null(postCdfTransform)) {
    postCdfTransform <- function(cdf) { 
      applyCdfGroups(cdf, lapply, function(group) {
        idxs <- group[[1]];
        if (is.null(dim(idxs))) {
          group;
        } else {
          list(indices=idxs[1,]);
        }
      })
    }
  }

  # Call main fit() function
  NextMethod("fit", this, ..., transform=transform, postCdfTransform=postCdfTransform);
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2006-09-10
# o Recreated.
############################################################################
