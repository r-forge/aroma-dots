# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE


.onAttach <- function(libname, pkgname) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Loading/installing affxparser
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Load 'affxparser'
  res <- suppressWarnings({
    require("affxparser");
  });

  # Not installed?
  if (!res) {
    if (interactive()) {
      cat("Package 'affxparser' is not available or could not be loaded. Will now try to install it from Bioconductor (requires working internet connection):\n");
      # To please R CMD check
      biocLite <- NULL; rm(biocLite);
      source("http://www.bioconductor.org/biocLite.R");
      biocLite("affxparser");
      # Assert that the package can be successfully loaded
      res <- require("affxparser");
      if (!res) {
        throw("Package 'affxparser' could not be loaded. Please install it from Bioconductor, cf. http://www.bioconductor.org/");
      }
    } else {
      warning("Package 'affxparser' could not be loaded. Without it ", pkgname, " will not work. Please install it from Bioconductor, cf. http://www.bioconductor.org/");
    }
  }
 

  pkg <- Package(pkgname);
  assign(pkgname, pkg, pos=getPosition(pkg));

  packageStartupMessage(getName(pkg), " v", getVersion(pkg), " (", 
     getDate(pkg), ") successfully loaded. See ?", pkgname, " for help.");
}
