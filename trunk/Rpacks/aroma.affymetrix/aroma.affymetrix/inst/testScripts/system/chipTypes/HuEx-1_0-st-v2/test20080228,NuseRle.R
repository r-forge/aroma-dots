library(aroma.affymetrix);

log <- Arguments$getVerbose(-3);
timestampOn(log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set (from previous file)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
name <- "Affymetrix-HeartBrain";
chipType <- "HuEx-1_0-st-v2";
cdfCore <- AffymetrixCdfFile$fromChipType(chipType, tags="core");
csTissue <- AffymetrixCelSet$fromName(name=name, cdf=cdfCore);
#Background correction and normalization
bc <- RmaBackgroundCorrection(csTissue);
csBCTissue <- process(bc);
qn <- QuantileNormalization(csBCTissue, typesToUpdate="pm");
csNTissue <- process(qn);
#Probe-level summarization
plmTissue <- ExonRmaPlm(csNTissue, mergeGroups=TRUE);
cesTissue <- getChipEffectSet(plmTissue);
plmNoMergeTissue <- ExonRmaPlm(csNTissue, mergeGroups=FALSE);
cesNoMergeTissue <- getChipEffectSet(plmNoMergeTissue);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Merged Groups:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Manually Calculate RLE/NUSE for first 100 units
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
units<-1:100
#Rle
theta <- extractMatrix(cesTissue,field="theta",units=units)
avg<-getAverageLog(cesTissue,field="intensities", mean="median",verbose=log)
thetaR <- extractMatrix(avg,field="theta",units=units)
RLE <- sweep(log2(theta),1,FUN="-",STATS=log2(thetaR))
#Nuse
theta <- extractMatrix(cesTissue, field="sdTheta",units=units)
avg<-getAverageLog(cesTissue,field="stdvs", mean="median", verbose=log)
thetaR <- extractMatrix(avg,field="theta",units=units)
NUSE <- sweep(log2(theta),1,FUN="/",STATS=log2(thetaR))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare to extractMatrix Results
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
testRLE<-extractMatrix(cesTissue,field="RLE",units=units)
stopifnot(identical(testRLE[,1], RLE[,1]));
testNUSE<-extractMatrix(cesTissue,field="NUSE",units=units)
stopifnot(identical(testNUSE[,1], NUSE[,1]));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare boxplot summaries of manual calculation to calculate{Rle|Nuse}BoxplotStats
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bxpRLE1<-boxplot.stats(RLE[,1])
testRLEBxp<-calculateRleBoxplotStats(cesTissue, arrays = 1:2, subset = units) 
stopifnot(identical(testRLEBxp[[1]]$stats, bxpRLE1$stats));

bxpNUSE1<-boxplot.stats(NUSE[,1])
testNUSEBxp<-calculateNuseBoxplotStats(cesTissue, arrays = 1:2, subset = units) 
stopifnot(identical(testNUSEBxp[[1]]$stats, bxpNUSE1$stats));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Repeat for No-Merged Groups:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Manually Calculate RLE/NUSE for first 100 units
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
units<-1:100
#Rle
theta <- extractMatrix(cesNoMergeTissue,field="theta",units=units)
avg<-getAverageLog(cesNoMergeTissue,field="intensities", mean="median",verbose=log)
thetaR <- extractMatrix(avg,field="theta",units=units)
RLE <- sweep(log2(theta),1,FUN="-",STATS=log2(thetaR))
#Nuse
theta <- extractMatrix(cesNoMergeTissue, field="sdTheta",units=units)
avg<-getAverageLog(cesNoMergeTissue,field="stdvs", mean="median", verbose=log)
thetaR <- extractMatrix(avg,field="theta",units=units)
NUSE <- sweep(log2(theta),1,FUN="/",STATS=log2(thetaR))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare to extractMatrix Results
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
testRLE<-extractMatrix(cesNoMergeTissue,field="RLE",units=units)
stopifnot(identical(testRLE[,1], RLE[,1]));
testNUSE<-extractMatrix(cesNoMergeTissue,field="NUSE",units=units)
stopifnot(identical(testNUSE[,1], NUSE[,1]));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Compare boxplot summaries of manual calculation to calculate{Rle|Nuse}BoxplotStats
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bxpRLE1<-boxplot.stats(RLE[,1])
testRLEBxp<-calculateRleBoxplotStats(cesNoMergeTissue, arrays = 1:2, subset = units) 
stopifnot(identical(testRLEBxp[[1]]$stats, bxpRLE1$stats));

bxpNUSE1<-boxplot.stats(NUSE[,1])
testNUSEBxp<-calculateNuseBoxplotStats(cesNoMergeTissue, arrays = 1:2, subset = units) 
stopifnot(identical(testNUSEBxp[[1]]$stats, bxpNUSE1$stats));
