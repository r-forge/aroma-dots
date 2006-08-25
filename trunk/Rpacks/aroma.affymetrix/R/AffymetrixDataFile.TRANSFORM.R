# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# 
#  P R O B E   S I G N A L   T R A N S F O R M A T I O N S
# 
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

setMethodS3("setTransformStatus", "AffymetrixDataFile", function(this, status=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  status <- Arguments$getLogical(status);


  oldStatus <- as.logical(this$.doTransform);
  this$.doTransform <- status;

  invisible(oldStatus);
}, protected=TRUE)

setMethodS3("transformOn", "AffymetrixDataFile", function(this, ...) {
  setTransformStatus(this, TRUE);
}, protected=TRUE)

setMethodS3("transformOff", "AffymetrixDataFile", function(this, ...) {
  setTransformStatus(this, FALSE);
}, protected=TRUE)

setMethodS3("isTransforming", "AffymetrixDataFile", function(this, ...) {
  as.logical(this$.doTransform);
}, protected=TRUE)

setMethodS3("hasTransforms", "AffymetrixDataFile", function(this, ...) {
  !is.null(this$.transforms);
}, protected=TRUE)

setMethodS3("setTransform", "AffymetrixDataFile", function(this, fcn, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Setting transformation function");
  this$.transforms <- list(fcn);
  verbose && exit(verbose);

  invisible(this);
}, protected=TRUE)


setMethodS3("addTransform", "AffymetrixDataFile", function(this, fcn, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Adding transformation function");

  transforms <- this$.transforms;
  transforms <- c(transforms, list(fcn));
  verbose && cat(verbose, "Number of functions afterwards: ", length(transforms));

  this$.transforms <- transforms;

  rm(transforms);

  verbose && exit(verbose);
}, protected=TRUE)




setMethodS3("transformProbeSignals", "AffymetrixDataFile", function(this, y, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Transforming probe signals");

  fcns <- this$.transforms;
  if (is.null(fcns))
    return(y);

  for (fcn in fcns) {
    y <- fcn(y, ...);
  }

  verbose && exit(verbose);

  y;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2006-05-15
# o Extracted from AffymetrixDataFile.R.
############################################################################
