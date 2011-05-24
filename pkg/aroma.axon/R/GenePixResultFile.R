setConstructorS3("GenePixResultFile", function(..., .verify=TRUE) {
  this <- extend(GenericDataFile(...), "GenePixResultFile",
    .fileHeader = NULL
  );

  if (.verify)
    verify(this, ..., verbose=verbose);
  this;
})


setMethodS3("as.character", "GenePixResultFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", this, ...);
  class <- class(s);
#  s <- c(s, sprintf("Number of text lines: %d", nbrOfLines(this, fast=TRUE)));

  class(s) <- class;
  s;
})

setMethodS3("getExtensionPattern", "GenePixResultFile", function(static, ...) {
  "[.](gpr)$";
}, static=TRUE)

setMethodS3("getPlatform", "GenePixResultFile", function(this, ...) {
  "Axon";
})

setMethodS3("getChipType", "GenePixResultFile", function(this, ...) {
  path <- getPath(this);
  dirname <- basename(path);
  dirname;
})


setMethodS3("verify", "GenePixResultFile", function(this, ..., verbose=FALSE) {
  # Nothing to do?
  if (is.null(getPathname(this)))
    return(invisible(this));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  return(this);

  verbose && enter(verbose, "Validating file contents");

  verbose && exit(verbose);

  invisible(this);
}, private=TRUE)



############################################################################
# HISTORY:
# 2010-05-23
# o Created.
############################################################################
