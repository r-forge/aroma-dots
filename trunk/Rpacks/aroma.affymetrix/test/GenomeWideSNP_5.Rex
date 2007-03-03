library(aroma.affymetrix)
source("init.R")

# How much more RAM too use than normal
ram <- 4

# Get data set
cs <- AffymetrixCelSet$fromName("Affymetrix_2007-GW_SNP5", chipType="GenomeWideSNP_5")

# Get the Affymetrix CSV genome information for this chip type
cdf <- getCdf(cs)
gi <- AffymetrixCsvGenomeInformation$fromChipType("GenomeWideSNP_5")
cdf$.gi <- gi

# Quantile normalize
qn <- QuantileNormalization(cs)
csN <- process(qn, verbose=verbose)

# Fit the PLM
plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE)
fit(plm, ram=ram, verbose=verbose)

# Get chip effects
ces <- getChipEffectSet(plm)

# Carry down the genome information
cdf <- getCdf(ces)
cdf$.gi <- gi

# Make sure average chip effects are calculated
getAverage(ces, force=TRUE, verbose=verbose)

# Setup the GLAD model
glad <- GladModel(ces)

# Let the ChromosomeExplorer create the results
ce <- ChromosomeExplorer(glad)
setArrays(ce, 1:5)
process(ce, chromosomes=9:6, verbose=verbose)
