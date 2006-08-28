savehistory();
closeAllConnections();

library(R.oo)
library(R.graphics)
library(R.utils)
library(affxparser)
library(aroma.apd)
library(aroma.affymetrix)

# Digest is still broken
#source("patches/digest.R")

# Patching during development
#source("~/braju.com.R/affxparser/affxparser/R/updateCel.R")
#source("~/braju.com.R/affxparser/affxparser/R/private.readCelHeaderV4.R")

tryCatch({
  opwd <- setwd("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/");
  source("rowSds.R")
  source("plotUtils.R")
  source("AffymetrixFile.R")
  source("AffymetrixFileSet.R")
  source("AffymetrixCdfFile.R")
  source("AffymetrixCelFile.R")
  source("AffymetrixCelFile.PLOT.R")
  source("AffymetrixCelSet.R")
  
  source("ParameterFile.R")
  source("ParameterCelFile.R")
  
  source("ProbeAffinityFile.R")
  source("ChipEffectFile.R")
  source("ChipEffectSet.R")
  source("AffymetrixUnitGroupsModel.R")
  source("ProbeLevelModel.R")
  
  source("LiWongProbeAffinityFile.R")
  source("AffymetrixLiWongModel.R")
  source("AffymetrixCnLiWongModel.R")
  
  source("RmaProbeAffinityFile.R")
  source("AffymetrixRmaModel.R")
  source("AffymetrixCnRmaModel.R")
  
  source("GdasAnnotationFile.R")
  source("GdasAnnotationSet.R")
}, finally = { setwd(opwd) })


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
