if (interactive()) savehistory();
library(aroma.affymetrix);
#source("FragmentLengthNormalization.R");
#source("UnitTypeScaleNormalization.R");
#source("DChipQuantileNormalization.R");
#source("AllelicCrosstalkCalibration.getSetsOfProbes.R");

verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full");
print(cdf);

csR <- AffymetrixCelSet$byName("HapMap270,6.0,CEU,founders", cdf=cdf);
print(csR);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Before starting, assert that all annotation files are available!
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
gi <- getGenomeInformation(cdf, type="UGP");
print(gi);

si <- getSnpInformation(cdf, type="UFL");
print(si);

#ufc <- AromaUfcFile$byChipType(getChipType(cdf));
#print(ufc);

acs <- AromaCellSequenceFile$byChipType(getChipType(cdf, fullname=FALSE));
print(acs);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calibrate for allelic crosstalk
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
version <- 1;
acc <- AllelicCrosstalkCalibration(csR, mergeShifts=TRUE, B=1, tags="*,v1b");
if (!isDone(acc)) {
  B <- getParameters(acc)$B;
  if (B <= 1) {
    sets1 <- getSetsOfProbes(acc, version=1, force=TRUE, verbose=verbose);
  }
  if (version == 4) {
    sets4 <- getSetsOfProbes(acc, version=4, force=TRUE, verbose=verbose);
    str(sets4);
    stopifnot(length(sets4$nonSNPs) > 0);
    if (B <= 1) stopifnot(identical(sets1$nonSNPs, sets4$nonSNPs));
  } else {
    if (B >= 1) stop("B >= 1 not supported: ", B);
    sets1 <- getSetsOfProbes(acc, version=1, force=TRUE, verbose=verbose);
    str(sets1);
    stopifnot(length(sets1$nonSNPs) > 0);
  }
}
print(acc);
csC <- process(acc, verbose=verbose);
print(csC);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalization for base-position effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bpn <- BasePositionNormalization(csC, target="zero", tags="*,z", shift=+300);
print(bpn);
csN <- process(bpn, verbose=verbose);
print(csN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Summary probe-level data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE);
print(plm);

# Fit CN probes quickly (~5-10s/array!)
if (length(findUnitsTodo(plm)) > 0) {
  units <- fitCnProbes(plm, verbose=verbose);
  str(units);
}

fit(plm, ram=5, verbose=verbose);
ces <- getChipEffectSet(plm);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize for PCR fragment-length effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
targets <- rep(list(function(...) log2(2200)), 2);
fln <- FragmentLengthNormalization(ces, targetFunctions=targets, tags="*,z");
print(fln);

cesN <- process(fln, verbose=verbose);
print(cesN);


# if (!getCombineAlleles(cesN)) {
#   as <- AlleleSummation(cesN);
#   print(as);
#   cesS <- process(as, verbose=verbose);
#   print(cesS);
#   cesN <- cesS;
# }

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize for (remaining) scale differences between unit types
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sn <- UnitTypeScaleNormalization(cesN, tags="*");
print(sn);
cesN2 <- process(sn, verbose=verbose);
print(cesN2);
