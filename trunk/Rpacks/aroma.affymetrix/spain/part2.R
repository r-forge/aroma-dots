##########################################################################
# Data set:
# GSE8605/
#   Mapping10K_Xba142/
#    GSM226867.CEL, ..., GSM226876.CEL [10 files]
# URL: http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE8605
##########################################################################
library("aroma.affymetrix");
log <- Arguments$getVerbose(-8, timestamp=TRUE);
par(pch=19); # Use large solid dots in all scatter plots

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chip type annotation data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType("Mapping50K_Hind240");
print(cdf);
getChipType(cdf);
getDimension(cdf);
nbrOfUnits(cdf);
nbrOfCells(cdf);

unitNames <- getUnitNames(cdf);
unit <- indexOf(cdf, "SNP_A-1745139");
unit;
getUnitNames(cdf, units=unit);


getPathname(cdf);
getFileSize(cdf);
getFilename(cdf);
getPath(cdf);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Genome build annotation data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
gi <- getGenomeInformation(cdf);
print(gi);

getData(gi, units=unit);



si <- getSnpInformation(cdf);
print(si);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sample data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
csR <- AffymetrixCelSet$byName("HapMap270,100K,CEU,3trios", cdf=cdf);
print(csR);

# Names, tags, full names
getName(csR);
getFullName(csR);
getTags(csR);

# Chip type
getChipType(csR);
getCdf(csR);

# Arrays
nbrOfArrays(csR);
getNames(csR);
getTimestamps(csR);

# Paths and filenames
getPath(csR);
getPathnames(csR);


# Look at array #3 in this data set
cf <- getFile(csR, 3);
print(cf);
getName(cf);
nbrOfCells(cf);

# All probes
y <- extractMatrix(cf, drop=TRUE);
str(y);

# PM probes
cellsPM <- getCellIndices(cdf, stratifyBy="pm", unlist=TRUE);
str(cellsPM);
cellsMM <- getCellIndices(cdf, stratifyBy="mm", unlist=TRUE);
str(cellsMM);

# Get the (PM,MM) signals
yPM <- extractMatrix(cf, cells=cellsPM, drop=TRUE);
yMM <- extractMatrix(cf, cells=cellsMM, drop=TRUE);

boxplot(data.frame(PM=yPM, MM=yMM), ylim=c(0,2^14));


# Look at the density of the PM signals for all arrays
plotDensity(csR, types="pm", ylim=c(0,0.35));

y <- extractMatrix(csR, cells=502:508);

# Extract all signals from ChrX
gi <- getGenomeInformation(cdf);
unitsOnChrX <- getUnitsOnChromosome(gi, 23);
cellsOnChrX <- getCellIndices(cdf, units=unitsOnChrX, unlist=TRUE);
y <- extractMatrix(csR, cells=cellsOnChrX);
boxplot(as.data.frame(log2(y)), ylim=c(0,16), las=2);

# QUESTIONS:
# Where are the raw data files stored?
# How many arrays/CEL files are there?

# Where is the CDF stored?


# Calculate residuals from PLM.
ae <- ArrayExplorer(csR);
setColorMaps(ae, "sqrt,yellow");
process(ae);
display(ae);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Crosstalk between similar probes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#acc <- AllelicCrosstalkCalibration(csR);
acc <- AllelicCrosstalkCalibration(csR, model="CRMAv2", tags="*,v2");
print(acc);
csC <- process(acc, verbose=log);
print(csC);

# Where are the calibrated CEL files stored?
# What new tags have been added to this calibrated data set?
# Now, quit R and redo the above.  What happens?

# With four different nucleotides, how many (heterozygote) allele pairs can there be?
# How many if you ignore the ordering?

# What does the PM signal densities look like?
plotDensity(csC, types="pm", ylim=c(0,0.35));

y <- extractMatrix(csC, cells=cellsOnChrX);
boxplot(as.data.frame(log2(y)), ylim=c(0,16), las=2);

array <- indexOf(csR, "NA06985");
plotAllelePairs(acc, array=array, pairs=1:6, what="input", xlim=c(0,2^15));
plotAllelePairs(acc, array=array, pairs=1:6, what="output", xlim=c(0,2^15));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sequence effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bpn <- BasePositionNormalization(csC, target="zero");
print(bpn);
csN <- process(bpn, verbose=log);
print(csN);

# What does the PM signal densities look like?
plotDensity(csN, types="pm", ylim=c(0,0.35));


# Allele-pairs
accDummy <- AllelicCrosstalkCalibration(csC);
array <- indexOf(csR, "NA06985");
plotAllelePairs(accDummy, array=array, pairs=1:6, what="input", xlim=c(0,2^15));



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe summarization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE);
print(plm);
fit(plm, verbose=log);

ces <- getChipEffectSet(plm);
print(ces);

# Where are the chip effect stored?
# Why are they called chip effects?


# Extracting chip effects (thetaA, thetaB)
cdf <- getCdf(ces);

unit <- indexOf(cdf, "SNP_A-1745139");
theta <- extractTheta(ces, units=unit, drop=TRUE);
rownames(theta) <- c("A", "B");
print(theta);

plot(log2(t(theta)), xlim=c(0,15), ylim=c(0,15));
abline(a=0, b=1);

freqB <- theta["B",]/(theta["A",]+theta["A",]);
genotypes <- rep("AB", length(freqB));
genotypes[freqB < 1/3] <- "AA";
genotypes[freqB > 2/3] <- "BB";

cols <- c(AA="red", BB="green", AB="blue");
plot(log2(t(theta)), col=cols[genotypes], pch=19, xlim=c(0,15), ylim=c(0,15));

# Easier way to get freqB estimates
data <- extractTotalAndFreqB(ces, units=unit, drop=TRUE);
freqB <- data["freqB",];

# Repeat for another SNP
unit <- indexOf(cdf, "SNP_A-1643368");
unit <- indexOf(cdf, "SNP_A-1652155");
theta <- extractTheta(ces, units=unit, drop=TRUE);
rownames(theta) <- c("A", "B");
freqB <- extractTotalAndFreqB(ces, units=unit, drop=TRUE)["freqB",];
genotypes <- rep("AB", length(freqB));
genotypes[freqB < 1/3] <- "AA";
genotypes[freqB > 2/3] <- "BB";
plot(log2(t(theta)), col=cols[genotypes], pch=19, xlim=c(0,15), ylim=c(0,15));
title(main=getUnitNames(cdf, units=unit));

ceR <- getAverageFile(ces);
thetaR <- extractTheta(ceR, units=unit, drop=TRUE);

#C <- 2*theta/thetaR;
#plot(t(C)), col=cols[genotypes], xlim=c(0,15), ylim=c(0,15));

qam <- QualityAssessmentModel(plm);
plotRle(qam);
plotNuse(qam);

n23 <- sapply(csR, getAttribute, "n23");
isXX <- (n23 == 2);
plotRle(qam, subset=unitsOnChrX);

# Calculate residuals from PLM.
rs <- calculateResidualSet(plm, verbose=log);
ae <- ArrayExplorer(rs);
setColorMaps(ae, c("log2,log2neg,rainbow", "log2,log2pos,rainbow"));
process(ae, interleaved="v");
display(ae);

# Do you find any spatial artifacts?
# Can you see them in the spatial plots of the raw intensities?


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PCR Fragment-length normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces, target="zero");
print(fln);
cesN <- process(fln, verbose=log);
print(cesN);

unit <- indexOf(cdf, "SNP_A-1643368");
unit <- indexOf(cdf, "SNP_A-1652155");
theta <- extractTheta(cesN, units=unit, drop=TRUE);
rownames(theta) <- c("A", "B");
freqB <- extractTotalAndFreqB(ces, units=unit, drop=TRUE)["freqB",];
genotypes <- rep("AB", length(freqB));
genotypes[freqB < 1/3] <- "AA";
genotypes[freqB > 2/3] <- "BB";
plot(log2(t(theta)), col=cols[genotypes], pch=19, xlim=c(0,15), ylim=c(0,15));
title(main=getUnitNames(cdf, units=unit));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# FreqB plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
array <- indexOf(cesN, "GSM226871");
ce <- getFile(cesN, array);
units <- getUnitsOnChromosome(gi, 2);
pos <- getPositions(gi, units);
freqB <- extractTotalAndFreqB(ce, units=units)[,"freqB"];
plot(pos, freqB, ylim=c(0,1), pch=19);
plot(density(freqB[pos < 90e6], adjust=0.3), col="blue", lwd=2);
lines(density(freqB[pos > 90e6], adjust=0.3), col="red", lwd=2);

# Clear deletion
array <- indexOf(cesN, "GSM226871");
ce <- getFile(cesN, array);
units <- getUnitsOnChromosome(gi, 9);
pos <- getPositions(gi, units);
freqB <- extractTotalAndFreqB(ce, units=units)[,"freqB"];
plot(pos, freqB, ylim=c(0,1), pch=19);
plot(density(freqB[pos < 40e6], adjust=0.3), col="blue", lwd=2);
lines(density(freqB[pos > 65e6], adjust=0.3), col="red", lwd=2);

# Clear deletion
array <- indexOf(cesN, "GSM226871");
ce <- getFile(cesN, array);
units <- getUnitsOnChromosome(gi, 13);
pos <- getPositions(gi, units);
freqB <- extractTotalAndFreqB(ce, units=units)[,"freqB"];
plot(pos, freqB, ylim=c(0,1), pch=19);
plot(density(freqB[pos < 40e6], adjust=0.3), col="blue", lwd=2);
lines(density(freqB[pos > 65e6], adjust=0.3), col="red", lwd=2);


array <- indexOf(cesN, "GSM226875");
ce <- getFile(cesN, array);
units <- getUnitsOnChromosome(gi, 5);
pos <- getPositions(gi, units);
freqB <- extractTotalAndFreqB(ce, units=units)[,"freqB"];
plot(pos, freqB, ylim=c(0,1), pch=19);
plot(density(freqB[pos < 45e6], adjust=0.3), col="blue", lwd=2);
lines(density(freqB[pos > 45e6], adjust=0.3), col="red", lwd=2);

Clear deletion
array <- indexOf(cesN, "GSM226876");
ce <- getFile(cesN, array);
units <- getUnitsOnChromosome(gi, 5);
pos <- getPositions(gi, units);
freqB <- extractTotalAndFreqB(ce, units=units)[,"freqB"];
plot(pos, freqB, ylim=c(0,1), pch=19);
plot(density(freqB[pos < 45e6], adjust=0.3), col="blue", lwd=2);
lines(density(freqB[pos > 45e6], adjust=0.3), col="red", lwd=2);


array <- indexOf(cesN, "GSM226874");
ce <- getFile(cesN, array);
units <- getUnitsOnChromosome(gi, 6);
pos <- getPositions(gi, units);
freqB <- extractTotalAndFreqB(ce, units=units)[,"freqB"];
plot(pos, freqB, ylim=c(0,1), pch=19);
plot(density(freqB[pos < 65e6], adjust=0.5), col="blue", lwd=2);
lines(density(freqB[pos > 65e6], adjust=0.5), col="red", lwd=2);


array <- indexOf(cesN, "GSM226874");
ce <- getFile(cesN, array);
units <- getUnitsOnChromosome(gi, 11);
pos <- getPositions(gi, units);
freqB <- extractTotalAndFreqB(ce, units=units)[,"freqB"];
plot(pos, freqB, ylim=c(0,1), pch=19);



ceR <- getAverage(cesN, verbose=log);





ce <- getFile(cesN, 1);
for (chr in getChromosomes(gi)) {
  units <- getUnitsOnChromosome(gi, chr);
  pos <- getPositions(gi, units=units) / 1e6; 
  thetaR <- extractTotalAndFreqB(ceR, units=units)[,"total"];
  data <- extractTotalAndFreqB(ce, units=units);
  data[,"total"] <- 2*data[,"total"] / thetaR;  

  fig <- sprintf("%s,Chr%d", getFullName(ce), chr);
  if (!devIsOpen(fig)) {
    devSet(fig);
    layout(matrix(1:2, ncol=1));
    par(mar=c(3,4,2,1)+0.1, pch=".");

    cn <- RawCopyNumbers(data[,"total"], pos, chromsome=chr);
    plot(cn, col="gray", cex=0.8, ylim=c(0,4));
    cnS <- gaussianSmoothing(cn, xOut=seq(xMin(cn), xMax(cn), by=1/2), sd=1);
    points(cnS, col="black");
    stext(side=3, pos=0, getName(ce));
    stext(side=3, pos=1, sprintf("Chr%d", chr));

    plot(pos, data[,"freqB"], cex=3, ylim=c(0,1));
    box(col="blue");
    stext(side=3, pos=0, getTags(cesN, collapse=","));

    devDone();
  }
}
