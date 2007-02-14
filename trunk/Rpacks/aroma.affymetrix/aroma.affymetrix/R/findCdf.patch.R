findCdf.patch <- function(chipType=NULL, paths=NULL, pattern="[.](c|C)(d|D)(f|F)$", ...) {
  settings <- getOption("affxparser.settings");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Customized search method?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  findFcn <- settings$methods$findCdf;
  if (!is.null(findFcn)) {
    pathnames <- findFcn(chipType=chipType, paths=paths, pattern=pattern, ...);
    if (!is.null(pathnames))
      return(pathnames);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup search path
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'paths':
  if (is.null(paths)) {
    paths <- paste(".", 
                   "cdf/", "data/cdf/", 
                   settings$paths$annotationData$cdf,
                   getOption("AFFX_CDF_PATH"), 
                   Sys.getenv("AFFX_CDF_PATH"),
             sep=";", collapse=";");
  }
  paths <- c(".", paths);
  paths <- unique(paths);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup search pattern
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(chipType)) {
    if (regexpr("[.](c|C)(d|D)(f|F)$", chipType) !=-1)
      warning("Argument 'chipType' of findCdf() has suffix '.cdf':", chipType);
    pattern <- paste(chipType, pattern, sep="");
  }

  findFiles(pattern=pattern, paths=paths, ...);
}

############################################################################
# HISTORY:
# 2007-02-12
# o Created.  Patched in patchAffxparser().
############################################################################  
