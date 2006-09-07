###########################################################################/**
# @RdocClass TotalCnAffinePlm
#
# @title "The TotalCnAffinePlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the affine model used by Bengtsson \& Hössjer (2006).
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
# @author
#
# \section{Model estimates}{
#   The estimated probe affinities are represented by the
#   @see "AffineProbeAffinityFile" class.  
# }
#
# \references{
# }
#*/###########################################################################
setConstructorS3("TotalCnAffinePlm", function(..., name="modelTotalCnAffine") {
  this <- extend(AffinePlm(..., name=name), "TotalCnAffinePlm")
  this <- setup(this);
  this;
})


setMethodS3("setup", "TotalCnAffinePlm", function(this, ...) {
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



setMethodS3("fit", "TotalCnAffinePlm", function(this, ..., transform=NULL, postCdfTransform=NULL) {
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
# 2006-08-28
# o Created from the corresponding RMA model.  More or less identical.
############################################################################
