#############################################################################
# Generic script to run all example files in the "inst" directory
#
# This script assumes the following directory structure:
#   <pkg>/
#    +-<pkg>/
#       +-incl/
#         +- foo.Rex (one example file)
#         +- bar.Rex (another example file)
#       +-man/
#       +-R/
#       +-DESCRIPTION
#    +-test/
#       +- runExamples.Rex  (this script)
#############################################################################
# Get current directory, i.e. <pkg>/test/
pwd <- getwd();

# Get the name of the parent directory, i.e. <pkg>
pkg <- basename(dirname(pwd));

# Assert that package path exists
pkgPath <- file.path("..", pkg)
if (!file.exists(pkgPath) || !file.info(pkgPath)$isdir)
  stop("Wrong directory structure");

# Assert that path exists
if (!all(c("DESCRIPTION", "R", "man", "incl") %in% list.files(pkgPath)))
  stop("Wrong directory structure");

# Load <pkg>
library(pkg, character.only=TRUE);

# List all example files in <pkg>/<pkg>/incl/
files <- list.files(pattern="[.](R|Rex)$", path=file.path(pkgPath, "incl"),
                                            all.files=TRUE, full.names=TRUE);

for (file in files) {
  cat("Running script '", basename(file), "'...\n", sep="");
  # Evaluated script in a anonymous function
  { function(file) {
      source(file, local=TRUE, echo=TRUE, max.deparse.length=2048);
  }}(file);
  gc();
  cat("Running script '", basename(file), "'...done\n", sep="");
}


#############################################################################
# HISTORY:
# 2006-03-30
# o Created.
#############################################################################
