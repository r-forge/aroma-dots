###########################################################################/**
# @RdocClass ParameterFile
#
# @title "The ParameterFile class"
#
# \description{
#  @classhierarchy
#
#  An ParameterFile object represents parameter estimates.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
# @visibility "private"
#*/###########################################################################
setConstructorS3("ParameterFile", function(..., parameters=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  extend(AffymetrixFile(...), "ParameterFile",
    parameters = parameters
  )
})

setMethodS3("nbrOfParameters", "ParameterFile", function(this, ...) {
  length(this$parameters);
})

setMethodS3("getParameters", "ParameterFile", function(this, ...) {
  this$parameters;
})

setMethodS3("getParameter", "ParameterFile", function(this, name, ...) {
  if (!name %in% this$parameters)
    throw("Unknown parameter: ", name);
  param <- this[[name]];
  param;
})



############################################################################
# HISTORY:
# 2006-05-31
# o Added abstract getExtension().
# 2006-04-01
# o Added getPath().
# 2006-03-17
# o Created.
############################################################################
