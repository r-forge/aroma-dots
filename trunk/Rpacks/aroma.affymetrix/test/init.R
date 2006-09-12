tryCatch({
  savehistory();
}, error = function(ex) {})

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

#source("../../affxparser/affxparser/R/convertCdf.R");
#source("../../affxparser/affxparser/R/writeCdf.R");
#source("../../affxparser/affxparser/R/applyCdfBlocks.R");
#sourceDirectory("../affxparser/affxparser/R/", recursive=FALSE);

sourceDirectory("../aroma.affymetrix/R/", recursive=FALSE);

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

