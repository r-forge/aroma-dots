setMethodS3("isFileModified", "default", function(filename, path=NULL, since=NULL, clear=FALSE, ...) {
  # Arguments 'filename' & 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path);

  timestamps <- getOption("fileModificationTimestamps");
  if (is.null(since)) {
    since <- timestamps[[pathname]];
  } else if (inherits(since, "POSIXct")) {
  } else {
    since <- as.POSIXct(since, ...);
  }

  # Check if file is modified later than 'since'
  timestamp <- file.info(pathname)$mtime;

  # Has been checked before?
  if (is.null(since) || clear) {
    isModified <- TRUE;
  } else {
    isModified <- (timestamp > since);
  }

  # Update known timestamps?
  if (isModified) {
    timestamps[[pathname]] <- timestamp;
    options(fileModificationTimestamps=timestamps);
  }

  isModified;
}) # isFileModified()


############################################################################
# HISTORY:
# 2008-07-24
# o Created.
############################################################################
