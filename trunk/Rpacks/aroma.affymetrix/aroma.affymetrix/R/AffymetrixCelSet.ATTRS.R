setMethodS3("setAttributesBy", "AffymetrixCelSet", function(this, object, ...) {
  methodName <- sprintf("setAttributesBy%s", class(object)[1]);
  if (!exists(methodName, mode="function")) {
    throw("No set function found: ", methodName);
  }
  
  fcn <- get(methodName, mode="function");
  tryCatch({
    fcn(this, object, ...);
  }, error = function(ex) {
    print(ex);
    throw("Failed to apply attributes by object of class: ", class(object)[1]);
  })
}, protected=TRUE)


setMethodS3("setAttributesBySampleAnnotationSet", "AffymetrixCelSet", function(this, sas, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  res <- lapply(sas, FUN=function(saf) {
    verbose && enter(verbose, "Applying sample annotations");
    on.exit({verbose && exit(verbose)});

    setAttributesBy(this, saf, ..., verbose=less(verbose));
  });
  invisible(res);
}, protected=TRUE)


setMethodS3("setAttributesBySampleAnnotationFile", "AffymetrixCelSet", function(this, saf, force=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setAttrs <- function(appliesTo, tags=NULL, ..., verbose=FALSE) { 
    verbose && enter(verbose, "Applying sample annotations");
    on.exit({verbose && exit(verbose)});

    args <- list(...);
    nargs <- length(args);
    if (!is.null(tags)) {
      tags <- tags[!is.na(tags)];
      tags <- tags[nchar(tags) > 0];
      if (length(tags) == 0)
        tags <- NULL;
    }

    if (!is.null(tags)) {
      # Split tags
      tags <- unlist(strsplit(tags, split=","), use.names=FALSE);
      tags <- trim(tags);
      verbose && cat(verbose, "Tags: ", paste(tags, collapse=", "));
      nargs <- nargs + 1;
    }

    # Nothing to do?
    if (nargs == 0)
      return();

    # Typically the below only applies to one sample
    verbose && cat(verbose, "Applies to ", length(appliesTo), " sample(s).");
    for (kk in seq(along=appliesTo)) { 
      idx <- appliesTo[kk];
      verbose && cat(verbose, "Sample: ", names(appliesTo)[kk]);

      # Get the CEL file
      cf <- getFile(this, idx);

      # Apply the attributes
      setAttributes(cf, ...);

      # Apply the tags
      if (!is.null(tags)) {
        setAttributesByTags(cf, tags);
      }
    };
  } # setAttrs()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'saf':
  if (!inherits(saf, "SampleAnnotationFile")) {
    throw("Argument 'saf' is not a SampleAnnotationFile: ", class(saf)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  names <- getFullNames(this);
  res <- apply(saf, names, FUN=setAttrs, force=force, verbose=verbose);

  invisible(res);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2007-03-14
# o Now setAttributesBySampleAnnotationFile() also set attributes.
# 2007-03-06
# o Added setAttributesBy().
# o Added setAttributesBySampleAnnotationFile().
# o Created.
############################################################################
