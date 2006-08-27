###########################################################################/**
# @RdocClass ProbeAffinityFile
#
# @title "The ProbeAffinityFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of probe affinities in probe-level models.
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




setMethodS3("getFirstCellIndices", "ProbeAffinityFile", function(this, units=NULL, ..., verbose=FALSE) {
  # Argument 'verbose': 
 verbose <- Arguments$getVerbose(verbose);

  res <- this$.firstCells;
  if (is.null(res)) {
    stratifyBy <- switch(this$model, pm="pm");
    cdf <- getCdf(this);
    res <- getFirstCellIndices(cdf, units=NULL, ..., stratifyBy=stratifyBy, verbose=verbose);
    this$.firstCells <- res;
  }

  # Subset?
  if (!is.null(units))
    res <- res[units];

  res;
}, protected=TRUE)


setMethodS3("findUnitsTodo", "ProbeAffinityFile", function(this, units=NULL, field="stdvs", ..., transform=NULL, verbose=FALSE) {
  # Argument 'field':
  field <- match.arg(field);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # Get the indices of the first cells in each unit group
  indices <- getFirstCellIndices(this, units=units, ...);

  if (!is.null(transform)) {
    indices <- transform(indices);
  }

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


setMethodS3("writeSpatial", "ProbeAffinityFile", function(this, ..., transform=NULL, zlim=c(0,3)) {
  if (is.null(transform)) {
    transform <- function(x) {
      # Probe-affinities are in (0,Inf]
      nok <- (x == 0);
      # Truncate zeros to smallest strictly positive value
      x[nok] <- min(x[!nok], na.rm=TRUE);

#      x <- log(x, base=2);
      x;
    }
  }

  NextMethod("writeSpatial", this, ..., transform=transform, zlim=zlim);
})

############################################################################
# HISTORY:
# 2006-08-26
# o Added writeSpatial().
# 2006-08-25
# o Added findUnitsTodo().
# o Added getFirstCellIndices(). Since reading all cell indices can take
#   a while it is cached in memory, but also on file (in case we restart).
# o Created from LiWongProbeAffinityFile.  The RMA version is almost 
#   identical so I made this a superclass of both.
############################################################################
