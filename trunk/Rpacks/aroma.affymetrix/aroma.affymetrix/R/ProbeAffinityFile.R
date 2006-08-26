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
#   An object of this class is typically obtain through the
#   \code{getProbeAffinities()} method for an 
#   @see "AffymetrixProbeAffinityModel" class.
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
    model = model
  )
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
  stratifyBy <- switch(this$model, pm="pm");
  NextMethod("updateUnits", this, data=data, ..., stratifyBy=stratifyBy);
}, protected=TRUE);


############################################################################
# HISTORY:
# 2006-08-25
# o Created from LiWongProbeAffinityFile.  The RMA version is almost 
#   identical so I made this a superclass of both.
############################################################################
