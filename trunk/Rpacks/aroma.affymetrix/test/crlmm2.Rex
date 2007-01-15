library(oligo);
library(aroma.affymetrix);
source("init.R");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
.Machine$float.eps <- sqrt(.Machine$double.eps);
tolerance <- .Machine$float.eps^0.5;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
path <- "rawData/Test100K/Mapping50K_Hind240";  # Works! 2006-12-13
path <- "rawData/Test500K/Mapping250K_Nsp";

path <- "rawData/SinclairA_etal_2006/Mapping250K_Sty";

ds <- AffymetrixCelSet$fromFiles(path);
print(ds);

cdf <- getCdf(ds);
pd <- PlatformDesign(cdf);

if (!exists("dsC")) {
  verbose && enter(verbose, "Calibration for PCR effects");
  sc <- PcrProbeCalibration(ds, flavor="oligo");
  dsC <- process(sc, verbose=verbose);
  verbose && exit(verbose);
  gc();
}

if (!exists("dsN")) {
  verbose && enter(verbose, "Quantile normalization");
  norm <- OligoQuantileNormalization(dsC);
  verbose && print(verbose, norm);
  dsN <- process(norm, verbose=verbose);
  verbose && print(verbose, dsN);
  verbose && exit(verbose);
  gc();
}


if (!exists("ces")) {
  verbose && enter(verbose, "Probe-level modelling");
  plm <- RmaSnpPlm(dsN, mergeStrands=FALSE, flavor="oligo");
  verbose && print(verbose, plm);
  unitsDone <- fit(plm, verbose=verbose);
  ces <- getChipEffects(plm);
  verbose && exit(verbose);
  gc();
}

if (!exists("calls", mode="list")) {
  verbose && enter(verbose, "Genotyping using CRLMM");
  crlmm <- Crlmm(ces);
  verbose && print(verbose, crlmm);
  calls <- fit(crlmm, verbose=verbose);
  verbose && str(verbose, calls);
  verbose && exit(verbose);
  gc();
}
