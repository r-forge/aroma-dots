savehistory();
closeAllConnections();
library(R.oo)
library(R.utils)
library(R.graphics)
library(R.cache)
library(affxparser)

#library(aroma.apd)
#library(aroma.affymetrix)

verbose <- Arguments$getVerbose(TRUE);

verbose && enter(verbose, "Sourcing all *.R files");

# Digest is still broken
#source("patches/digest.R")

# Patching during development
#source("~/braju.com.R/affxparser/affxparser/R/updateCel.R")
#source("~/braju.com.R/affxparser/affxparser/R/private.readCelHeaderV4.R")

tryCatch({
  sourceDirectory("../aroma.affymetrix/R/", recursive=FALSE);
}, finally = { setwd(opwd) })

verbose && exit(verbose);

# Setup up the search path to ImageMagick
imageMagickConvert <- function(srcfile, destfile, format, options=NULL, ...) {
  pathname <- "C:/Program Files/ImageMagick-6.2.7-Q16/convert";
  # Escape everything
  pathname <- sprintf("\"%s\"", pathname);
  srcfile <- sprintf("\"%s\"", srcfile);
  destfile <- sprintf("\"%s\"", destfile);
  system(paste(pathname, options, srcfile, destfile));
}
options(imageConverter=imageMagickConvert);

