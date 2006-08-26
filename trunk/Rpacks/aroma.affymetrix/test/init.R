savehistory();
closeAllConnections();

library(R.oo)
library(R.utils)
library(affxparser)
library(aroma.apd)
library(aroma.affymetrix)

# Patching during development
#source("~/braju.com.R/affxparser/affxparser/R/updateCel.R")
#source("~/braju.com.R/affxparser/affxparser/R/private.readCelHeaderV4.R")

source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixFileSet.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixCdfFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixCelFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixCelSet.R")

source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/ParameterFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/ParameterCelFile.R")

source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/ProbeAffinityFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixUnitGroupsModel.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/ProbeLevelModel.R")

source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/LiWongProbeAffinityFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixLiWongModel.R")

source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/RmaProbeAffinityFile.R")
source("~/braju.com.R/aroma.affymetrix/aroma.affymetrix/R/AffymetrixRmaModel.R")
