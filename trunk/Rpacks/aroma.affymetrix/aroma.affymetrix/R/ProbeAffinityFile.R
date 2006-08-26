###########################################################################/**
# @RdocClass ProbeAffinityFile
#
# @title "The ProbeAffinityFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of probe affinities in the 
#  RMA model.
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
# @author
# 
# \seealso{
#   An object of this class is typically obtained through the
#   \code{getProbeAffinities()} method for the @see "ProbeLevelModel" class.
# }
#
#*/###########################################################################
setConstructorS3("ProbeAffinityFile", function(..., model=c("pm")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  extend(ParameterCelFile(...), "ProbeAffinityFile",
    "cached:.firstCells" = NULL,
    model = model
  )
})


# Note: This method 
setMethodS3("getFirstCellIndices", "ProbeAffinityFile", function(this, units=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  res <- this$.firstCells;
  if (is.null(res)) {
    cdf <- getCdf(this);
    stratifyBy <- switch(this$model, pm="pm");
    verbose && enter(verbose, "Checking in file cache");
    key <- list(chipType=getChipType(cdf), stratifyBy=stratifyBy);
    res <- loadCache(key=key);
    verbose && exit(verbose);
    if (is.null(res)) {
      verbose && enter(verbose, "Reading all cell indices");
      res <- getCellIndices(cdf, units=NULL, ..., stratifyBy=stratifyBy);
      verbose && exit(verbose);
    
      verbose && enter(verbose, "Extracting the first cell in each unit group");
      # For each unit and each group, get the index of the first cell.
      res <- applyCdfGroups(res, function(groups) {
        # For each group, pull out the first cell.
        lapply(groups, FUN=function(group) {
          # group$indices[1] == group[[1]][1] == ...
          list(indices=.subset(.subset2(group, 1), 1));
        })
      });
      verbose && exit(verbose);
      saveCache(key=key, res);
    }

    # Storing in cache
    this$.firstCells <- res;
  }

  # Subset?
  if (!is.null(units)) {
    res <- res[units];
  }

  res;
}, protected=TRUE)


setMethodS3("findUnitsTodo", "ProbeAffinityFile", function(this, units=NULL, field="stdvs", ..., verbose=FALSE) {
  # Argument 'field':
  field <- match.arg(field);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # Get the indices of the first cells in each unit group
  indices <- getFirstCellIndices(this, units=units, ...);

  # Keep only the first group
  indices <- applyCdfGroups(indices, .subset, 1);
  
  # Flatten to get a vector
  indices <- unlist(indices, use.names=FALSE);

  # Get the 'stdvs' for these cells
  readIntensities <- FALSE;
  readStdvs <- FALSE;
  readPixels <- FALSE;
  if (field == "stdvs") {
    readStdvs <- TRUE;
  }
  value <- readCel(getPathname(this), indices=indices, readIntensities=readIntensities, readStdvs=readStdvs, readPixels=readPixels);
  value <- value[[field]];

  # Identify units for which the stdvs <= 0.
  todo <- which(value <= 0);
  if (!is.null(units))
    todo <- units[todo];

  todo;
})


setMethodS3("createFrom", "ProbeAffinityFile", function(static, ..., filename="probeAffinities.CEL") {
  createFrom.ParameterCelFile(static, ..., filename=filename);
}, static=TRUE);



setMethodS3("readUnits", "ProbeAffinityFile", function(this, ...) {
  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  stratifyBy <- switch(this$model, pm="pm");
  NextMethod("readUnits", this, ..., stratifyBy=stratifyBy);
});


setMethodS3("updateUnits", "ProbeAffinityFile", function(this, data, ...) {
  # Note that the actually call to the encoding is done in updateUnits()
  # of the superclass.
  NextMethod("updateUnits", this, data=data, ...);
}, protected=TRUE);


############################################################################
# HISTORY:
# 2006-08-25
# o Added findUnitsTodo().
# o Added getFirstCellIndices(). Since reading all cell indices can take
#   a while it is cached in memory, but also on file (in case we restart).
# o Created from LiWongProbeAffinityFile.  The RMA version is almost 
#   identical so I made this a superclass of both.
############################################################################
