##########################################################################
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup annotation data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType("Mapping50K_Hind240");
print(cdf);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup raw data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName("HapMap270,100K,CEU,3trios", cdf=cdf);
print(csR);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Correct for crosstalk between alleles
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2", tags="*,v2");
print(acc);
csC <- process(acc, verbose=log);
print(csC);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Correct for probe-sequence dependent effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bpn <- BasePositionNormalization(csC, target="zero");
print(bpn);
csN <- process(bpn, verbose=log);
print(csN);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Summarize probesets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE);
print(plm);
fit(plm, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Correct for PCR fragment-length effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ces <- getChipEffectSet(plm);
print(ces);
fln <- FragmentLengthNormalization(ces, target="zero");
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get total chip effect estimates, i.e. theta = thetaA + thetaB
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cbs <- CbsModel(cesN);
print(cbs);

ce <- ChromosomeExplorer(cbs);
process(ce, verbose=log);