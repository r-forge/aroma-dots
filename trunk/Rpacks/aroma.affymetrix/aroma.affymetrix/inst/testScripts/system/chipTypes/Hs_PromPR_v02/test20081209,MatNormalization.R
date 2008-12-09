library("aroma.affymetrix");
## source("../aroma.affymetrix/R/AffymetrixCelSet.convertToUnique.R");
## source("../aroma.affymetrix/R/AromaCellPositionFile.R");
log <- Arguments$getVerbose(-20, timestamp=TRUE);

cdf <- AffymetrixCdfFile$fromChipType("Hs_PromPR_v02", tags="Harvard,ROIs");
print(cdf);

csR <- AffymetrixCelSet$fromName("MNtest", cdf=cdf);
print(csR);

mn <- MatNormalization(csR, numChunks=20);
f <- process(mn, verbose=-50);

# Get loci on chromosome 1
u <- which( substr(getUnitNames(cdf),1,5) == "chr1F" );
ind <- getCellIndices(cdf, units=u, stratifyBy="pm", unlist=TRUE, useNames=FALSE);
d <- extractMatrix(f, cells=ind, verbose=log);

# test the conversion to the unique CDF
cdfU <- getUniqueCdf(cdf, verbose=more(log, 60));
csU <- convertToUnique(f, verbose=log);

# Get loci on chromosome 1
u <- which( substr(getUnitNames(cdfU),1,5) == "chr1F" );
indU <- getCellIndices(cdfU, units=u, stratifyBy="pm", unlist=TRUE, useNames=FALSE);
dU <- extractMatrix(csU, cells=indU, verbose=log);

stopifnot(all.equal(dU[,1], d[,1]));

acp <- AromaCellPositionFile$byChipType(getChipType(cdfU))

ch <- acp[indU,][,1];
sp <- acp[indU,][,2];

keyAroma <- paste("chr", ch, ".", sp, sep="");

barTxtFile <- read.table("TESTIP.bar.txt", sep="\t", header=FALSE, comment.char="#", nrow=450000);
keyBar <- paste(barTxtFile$V1, barTxtFile$V2, sep=".");

common <- intersect(keyAroma, keyBar);
stopifnot(length(common) == nrow(d));

mA <- match(common, keyAroma);
mB <- match(common, keyBar);

stopifnot( mean( (log2(d[mA,1])-barTxtFile[mB,3])^2 ) < .001 );
