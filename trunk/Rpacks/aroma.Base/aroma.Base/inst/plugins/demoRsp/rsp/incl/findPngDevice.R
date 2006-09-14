####################################################################
# The findPngDevice() returns a device that generated PNG files,
# or NULL if no such device exists.  This function is used because
# png() is not available on Unix if R runs in batch mode (X11 is 
# missing), and therefore bitmap() is returned, if Ghostscript is
# available.
####################################################################
findPngDevice <- function(...) {
  file <- tempfile("hasPng");
  on.exit({
    if (file.exists(file))
      file.remove(file);
  })

  # First, try png() device, if it exists...
  if (capabilities("png")) {
    tryCatch({
      png(file);
      plot(0);
      dev.off();
      if (file.exists(file))
        return(png);
    }, error = function(error) {
    });
  }

  # Second, try bitmap() device, which utilizes ghostscript...
  tryCatch({
    bitmap(file);
    plot(0);
    dev.off();
    if (file.exists(file))
      return(bitmap);
  }, error = function(error) {
  });

  # Third, try bitmap() device again, but search for ghostview first
  # (if found, environment variable R_GSCMD is updated too).
  gscmd <- System$findGhostscript();
  if (is.null(gscmd))
    return(NULL);

  tryCatch({
    bitmap(file);
    plot(0);
    dev.off();
    if (file.exists(file))
      return(bitmap);
  }, error = function(error) {
  });
  
  NULL;
}
