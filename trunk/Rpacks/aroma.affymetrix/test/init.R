savehistory();
closeAllConnections();

library(R.oo)
library(R.utils)
library(affxparser)
library(aroma.apd)
library(aroma.affymetrix)

# Digest is still broken
#source("patches/digest.R")

source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/rowSds.R")

# Patching during development
#source("~/braju.com.R/affxparser/affxparser/R/updateCel.R")
#source("~/braju.com.R/affxparser/affxparser/R/private.readCelHeaderV4.R")

source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixFileSet.R")
# source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixFileSet.SUBSETTING.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixCdfFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixCelFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixCelFile.PLOT.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixCelSet.R")

source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/ParameterFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/ParameterCelFile.R")

source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/ProbeAffinityFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/ChipEffectFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/ChipEffectSet.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixUnitGroupsModel.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/ProbeLevelModel.R")

source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/LiWongProbeAffinityFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixLiWongModel.R")

source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/RmaProbeAffinityFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixRmaModel.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixCnRmaModel.R")

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
