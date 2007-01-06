###########################################################################/**
# @RdocClass AffymetrixPowerTools
#
# @title "The AffymetrixPowerTools class"
#
# \description{
#  @classhierarchy
#
#  This class represents an interface to the Affymetrix Power Tools (APT)
#  applications [1].
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#  Note that the methods of this class are calling the APT [1] binaries,
#  which must be installed on the system.
# }
#
# @author
#
# \seealso{
# }
#
# \references{
#  [1] Affymetrix, Affymetrix Power Tools (APT) software, Dec 2006.
#      \url{http://www.affymetrix.com/support/developer/powertools/index.affx}
#      \cr
# }
#
# @keyword internal
#*/###########################################################################
setConstructorS3("AffymetrixPowerTools", function(...) {
  extend(Object(), "AffymetrixPowerTools",
    ...
  )
}, private=TRUE)


setMethodS3("writeChrXFile", "AffymetrixPowerTools", function(static, cdf, path=NULL, overwrite=FALSE, ...) {
  # Argument 'cdf':
  if (!inherits(cdf, "AffymetrixCdfFile")) {
    throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  # Argument 'path':
  path <- Arguments$getWritablePath(path);

  # Retrieve the genome information
  gi <- getGenomeInformation(cdf);

  # Check if output file can be generated
  chipType <- getChipType(cdf);
  filename <- sprintf("%s.chrx", chipType);
  pathname <- filePath(path, filename);
  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=!overwrite);

  # Get the names of all units on chromosome X
  units <- getUnitsOnChromosome(gi, chromosome="X");
  units <- getUnitNames(cdf, units=units);

  # Write to file  
  lines <- c("all_chrx_no_par", units);
  writeLines(con=pathname, lines);

  invisible(pathname);
}, static=TRUE)


setMethodS3("as.character", "AffymetrixPowerTools", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)





############################################################################
# HISTORY:
# 2006-12-01
# o Added writeChrXFile().
# o Created a skeleton.
############################################################################
