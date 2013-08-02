# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE


.onAttach <- function(libname, pkgname) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  retrieveAffxparser <- function(pkgname, ...) {
    # Trick 'R CMD check' to not generate NOTEs.
    requireX <- base::require;
    catX <- base::cat;

    reqPkgName <- "affxparser";
    res <- suppressWarnings({
      requireX(reqPkgName, character.only=TRUE);
    });

    # Not installed?
    if (!res) {
      if (interactive()) {
        # Trick 'R CMD check' to not generate NOTEs.
        catX("Package 'affxparser' is not available or could not be loaded. Will now try to install it from Bioconductor (requires working internet connection):\n");
        # To please R CMD check
        biocLite <- NULL; rm(list="biocLite");
        source("http://www.bioconductor.org/biocLite.R");
        biocLite(reqPkgName);
        # Assert that the package can be successfully loaded
        res <- requireX(reqPkgName, character.only=TRUE);
        if (!res) {
          throw("Package 'affxparser' could not be loaded. Please install it from Bioconductor, cf. http://www.bioconductor.org/");
        }
      } else {
        warning("Package 'affxparser' could not be loaded. Without it ", pkgname, " will not work. Please install it from Bioconductor, cf. http://www.bioconductor.org/");
      }
    }
  } # retrieveAffxparser()


  # Loading/installing affxparser
  retrieveAffxparser(pkgname);

  pkg <- Package(pkgname);
  assign(pkgname, pkg, pos=getPosition(pkg));

  startupMessage(pkg);
}


############################################################################
# HISTORY:
# 2012-04-01
# o Now .onAttach() hides the fact that it tries to load 'affxparser'
#   from 'R CMD check'.
############################################################################
