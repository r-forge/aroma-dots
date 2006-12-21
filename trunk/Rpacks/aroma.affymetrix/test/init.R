# Start R with nohup on shadowfax:
#  echo 'source("models.R")' | nohup R-2.4.0 --save &

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup up search paths
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
options("APT_PATH"="/home/users/lab0605/bengtsson/shared/Bioinformatics/software/AffymetrixPowerTools/bin/linux_x86_64");

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

svnCommit <- function(...) {
  if (atHome)
    return(FALSE);
  opwd <- getwd();
  on.exit(setwd(opwd));
  system("svn commit -m \"\"");
  setwd("../aroma.affymetrix/");
  system("svn commit -m \"\"");
  TRUE;
} # svnCommit()

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


glog2 <- function(x, scale=log2(2^16)/asinh(2^16), ...) {
  scale*(asinh(x)-asinh(1));
} # glog2()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tryCatch({
  savehistory();
}, error = function(ex) {})

closeAllConnections();


library(R.oo)
library(R.utils)

verbose <- Arguments$getVerbose(-8);
timestampOn(verbose);

library(R.graphics)
library(R.cache)
library(aroma.apd)
library(affxparser)
library(digest)
library(geneplotter)
library(RColorBrewer)
library(R.rsp)

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
} else {
  options("R.cache.path"="/tmp/hb/.Rcache/");
  mkdirs(getOption("R.cache.path"));
}

sourceDirectory("patches/", recursive=FALSE, modifiedOnly=TRUE);
sourceDirectory("../aroma.affymetrix/R/", recursive=FALSE, modifiedOnly=TRUE);

verbose && exit(verbose);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Change to a nicer palette
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cols <- brewer.pal(8, "Dark2");
palette(cols);

NORMAL <- 1; 
BOLD <- 2; 
ITALIC <- 3; 
BOLDITALIC <- 4;
