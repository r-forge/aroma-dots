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
#   \item{probeModel}{The specific type of probe model.}
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
setConstructorS3("ProbeAffinityFile", function(..., probeModel=c("pm", "mm", "pm-mm", "min1(pm-mm)")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'probeModel':
  probeModel <- match.arg(probeModel);

  extend(ParameterCelFile(...), "ProbeAffinityFile",
    "cached:.firstCells" = NULL,
    probeModel = probeModel
  )
})

setMethodS3("clearCache", "ProbeAffinityFile", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c(".firstCells")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)


setMethodS3("getCellIndices", "ProbeAffinityFile", function(this, ...) {
  stratifyBy <- switch(this$probeModel, "pm"="pm", "mm"="mm", "pm-mm"="pm", "min1(pm-mm)"="pm");
  cdf <- getCdf(this);
  getCellIndices(cdf, ..., stratifyBy=stratifyBy);
})


setMethodS3("readUnits", "ProbeAffinityFile", function(this, units=NULL, cdf=NULL, ...) {
  if (is.null(cdf))
    cdf <- getCellIndices(this, units=units);

  # Note that the actually call to the decoding is done in readUnits()
  # of the superclass.
  NextMethod("readUnits", this, cdf=cdf, readStdvs=TRUE, readPixels=TRUE, ...);
});


setMethodS3("updateUnits", "ProbeAffinityFile", function(this, units=NULL, cdf=NULL, data, ...) {
  if (is.null(cdf))
    cdf <- getCellIndices(this, units=units);

  # Note that the actually call to the encoding is done in updateUnits()
  # of the superclass.
  NextMethod("updateUnits", this, cdf=cdf, data=data, ...);
}, private=TRUE)




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
}, private=TRUE)



############################################################################
# HISTORY:
# 2007-01-03
# o Renamed constructor argument 'model' to 'probeModel'.
# 2006-09-11
# o Update read- and updateUnits() to make use of getCellIndices().
# o Added getCellIndices().
# 2006-08-26
# o Added writeSpatial().
# 2006-08-25
# o Added findUnitsTodo().
# o Added getFirstCellIndices(). Since reading all cell indices can take
#   a while it is cached in memory, but also on file (in case we restart).
# o Created from LiWongProbeAffinityFile.  The RMA version is almost 
#   identical so I made this a superclass of both.
############################################################################
