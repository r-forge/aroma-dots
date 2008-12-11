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

samples <- c("NA06993", "NA06985", "NA06991", "NA06994", "NA07000", "NA07029", "NA12005", "NA12006", "NA10839");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup annotation data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType("Mapping50K_Hind240");
print(cdf);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup raw data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName("HapMap270,100K,CEU,3trios", cdf=cdf);
# Interpret '_' (underscores) as tag separators (',')
setFullNamesTranslator(csR, function(names) gsub("_", ",", names));
#keep <- indexOf(csR, samples);
#csR <- extract(csR, keep);
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
plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE);
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
as <- AlleleSummation(cesN);
print(as);
cesT <- process(as, verbose=log);
print(cesT);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Find CN regions using the CBS segmentation method
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cbs <- CbsModel(cesT);
print(cbs);

ce <- ChromosomeExplorer(cbs);
process(ce, arrays=1, chromosomes=18, verbose=log);
process(ce, verbose=log);


# Look at:
# Chr23: All samples
# Chr22: NA06985, NA06991, NA06994, NA07000
