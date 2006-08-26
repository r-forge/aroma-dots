setMethodS3("transformAffine", "AffymetrixCelFile", function(this, outPath=file.path("transAffine", getChipType(this)), offset=0, scale=1, subsetToUpdate=NULL, typesToUpdate=NULL, ..., overwrite=FALSE, skip=!overwrite, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'outPath':
  outPath <- Arguments$getWritablePathname(outPath);
  if (identical(getPath(this), outPath)) {
    throw("Cannot not transform data. Argument 'outPath' refers to the same path as the path of the data file to be transformed: ", outPath);
  }

  # Argument 'offset':
  offset <- Arguments$getDouble(offset);

  # Argument 'scale':
  scale <- Arguments$getDouble(scale, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generating output pathname
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- basename(getPathname(this));
  pathname <- Arguments$getWritablePathname(filename, path=outPath, 
                                         mustNotExist=(!overwrite && !skip));

  # Already shifted?
  if (isFile(pathname) && skip) {
    verbose && cat(verbose, "Transformed data file already exists: ", pathname);
    return(fromFile(this, pathname));
  }

  # Get probe signals
  x <- getFields(this, fields="intensities", ..., verbose=verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);
  subsetToUpdate <- identifyCells(cdf, probes=subsetToUpdate, types=typesToUpdate, verbose=verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Shift intensities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, sprintf("Transforming probe intensities by (offset,scale)=(%.1f,%.2f) ", offset, scale));
  x[subsetToUpdate] <- offset + scale*x[subsetToUpdate];
  rm(subsetToUpdate);
  verbose && exit(verbose);

  # Write normalized data to file
  verbose && enter(verbose, "Writing transformed probe signals");

  # Copy CEL file and update the copy
  copyCel(from=getPathname(this), to=pathname, overwrite=overwrite);
  updateCel(pathname, intensities=x);

  verbose && exit(verbose);

  # Return transformed data file object
  fromFile(this, pathname);
}) # transformAffine()



############################################################################
# HISTORY:
# 2006-08-25
# o Move to class AffymetrixCelFile and output is now CEL files only.
# 2006-07-28
# o Added transformAffine().
############################################################################
