# Start R with nohup on shadowfax:
#  echo 'source("models.R")' | nohup R-2.4.0 --save &

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
atHome <- (Sys.getenv("HOST") == "");

svnUpdate <- function(...) {
  if (atHome)
    return(FALSE);
  opwd <- getwd();
  on.exit(setwd(opwd));
  system("svn update");
  setwd("../aroma.affymetrix/");
  system("svn update");
  TRUE;
} # svnUpdate()

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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tryCatch({
  savehistory();
}, error = function(ex) {})

closeAllConnections();


library(R.oo)
library(R.utils)

verbose <- Arguments$getVerbose(-3);
timestampOn(verbose);

library(R.graphics)
library(R.cache)
library(aroma.apd)
library(affxparser)
library(digest)
library(geneplotter)
library(RColorBrewer)

#library(aroma.apd)
#library(aroma.affymetrix)

if (!atHome)
  svnUpdate();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Patch/source
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
verbose && enter(verbose, "Sourcing all *.R files");

if (atHome) {
  source("~/braju.com.R/R.utils/R.utils/R/sourceDirectory.R");
  source("../../affxparser/affxparser/R/findCdf.R");
}
source("Verbose.R");

#source("../../affxparser/affxparser/R/convertCdf.R");
#source("../../affxparser/affxparser/R/writeCdf.R");
#source("../../affxparser/affxparser/R/applyCdfBlocks.R");
#sourceDirectory("../affxparser/affxparser/R/", recursive=FALSE);

sourceDirectory("../aroma.affymetrix/R/", recursive=FALSE);

verbose && exit(verbose);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Change to a nicer palette
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cols <- brewer.pal(8, "Dark2");
palette(cols);

