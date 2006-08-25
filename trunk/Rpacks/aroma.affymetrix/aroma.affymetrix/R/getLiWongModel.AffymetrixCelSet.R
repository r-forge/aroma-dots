setMethodS3("getLiWongModel", "AffymetrixCelSet", function(this, path="modelLiWong", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'path':
  path <- Arguments$getWritablePath(path=path);

  filename <- "probeAffinities.cel";
  pathname <- Arguments$getWritablePathname(filename, path=path);
  if (!isFile(pathname)) {
    verbose && enter(verbose, "Creating an empty CEL file for probe-affinity estimates.");
    file.copy(getPathname(ds[[1]]), pathname);
    # Clear the file
    updateCel(pathname, pixels=rep(0:0, length=nbrOfCells(getCdf(this))));
    verbose && exit(verbose);
  }

  AffymetrixLiWongModel(path=path, ..., dataset=this);
})
