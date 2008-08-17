if (interactive()) savehistory();
library(aroma.affymetrix);
#source("FragmentLengthNormalization.R");
source("UnitTypeScaleNormalization.R");
#source("DChipQuantileNormalization.R");

verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cdf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6", tags="Full");
print(cdf);

cs <- AffymetrixCelSet$byName("HapMap270,6.0,CEU,founders", cdf=cdf);
print(cs);


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
acc <- AllelicCrosstalkCalibration(cs);
print(acc);
csC <- process(acc, verbose=verbose);
print(csC);

## bc <- RmaBackgroundCorrection(cs);
## print(bc);
## csBC <- process(bc, verbose=verbose);
## print(csBC);
## csC <- csBC;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalization for base-count effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# bcn <- BaseCountNormalization(csC, shift=+300);
# print(bcn);
# csN <- process(bcn, verbose=verbose);
# print(csN);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalization for base-position effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bpn <- BasePositionNormalization(csC, target="zero", tags="*,z", shift=+300);
bpn <- BasePositionNormalization(cs, target="zero", tags="*,z");
print(bpn);
csN <- process(bpn, verbose=verbose);
print(csN);

acc <- AllelicCrosstalkCalibration(csN);
print(acc);
csC <- process(acc, verbose=verbose);
print(csC);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Summary probe-level data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE);
plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=TRUE);
plm <- RmaCnPlm(csC, mergeStrands=TRUE, combineAlleles=TRUE, shift=+300);
#plm <- AvgCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE);
plm <- AvgCnPlm(csC, mergeStrands=TRUE, combineAlleles=FALSE);
print(plm);

fit(plm, ram=5, verbose=verbose);
ces <- getChipEffectSet(plm);
print(ces);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize for PCR fragment-length effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fln <- FragmentLengthNormalization(ces, shift=300);
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
  
#cesN <- cesN2;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize for fragment equivalent class effects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#fec <- FragmentEquivalentClassNormalization(cesN);
#print(fec);

#cesN2 <- process(fec, verbose=verbose);  # ~30s/array
#print(cesN2);

##sn <- ScaleNormalization(cesN);
##cesN2 <- process(sn, verbose=verbose);
##print(cesN2);
##cesN <- cesN2;

#if (!exists("cesN")) {
#  cdfM <- getMonoCell(cdf);
#  cesN <- CnChipEffectSet$byName("HapMap270,6.0,CEU,founders", tags="ACC,ra,-XY,AVG,+300,A+B,FLN,-XY", cdf=cdfM, mergeStrands=TRUE, combineAlleles=TRUE)
#}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get estimates on ChrX
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- sprintf("%s,units23.RData", getChipType(gi));
if (isFile(pathname)) {
  units23 <- loadObject(pathname);
} else {
  units23 <- getUnitsOnChromosome(gi, chromosome=23);
  saveObject(units23, file=pathname);
}
str(units23);  # 87204 units

path <- "rocData";
mkdirs(path);
pathname <- filePath(path, sprintf("%s,units23.RData", getFullName(cesN)));
if (isFile(pathname)) {
  theta <- loadObject(pathname);
} else {
  data <- extractTotalAndFreqB(cesN, units=units23, drop=TRUE);
  str(data);
  theta <- data[,1,];
  str(theta);
  saveObject(theta, file=pathname);
}
print(getFullName(cesN));

# Write (total, freqB) files to cnData/, if missing
setList <- getTotalAndFreqBSets(cesN, verbose=verbose);
print(setList);

path <- "chrData";
path <- file.path(path, getFullName(cesN), getChipType(cesN, fullname=FALSE));
mkdirs(path);
filename <- sprintf("%s,units23,totalFreqB.RData", getFullName(cesN));
pathname <- filePath(path, filename);
if (isFile(pathname)) {
  theta <- loadObject(pathname);
} else {
  arrayNames <- getNames(setList[[1]]);
  dimnames <- list(NULL,  dim=names(setList), arrayNames);
  dim <- sapply(dimnames, FUN=length);
  dim[1] <- length(units23);
  data <- array(as.double(NA), dim=dim, dimnames=dimnames);
  for (what in names(sets)) {
    data[,what,] <- extractMatrix(sets[[what]], rows=units23);
  }
  str(data);
  saveObject(data, file=pathname);
}
print(getFullName(cesN));

stop();

rownames(theta) <- getUnitNames(cdf, units=units23);
str(theta);

thetaR <- rowMedians(theta, na.rm=TRUE);
M <- log2(theta/thetaR);


isSnp <- (regexpr("SNP_", rownames(theta)) != -1); # 38004 SNPs
isCnp <- (regexpr("CN_", rownames(theta)) != -1);  # 49200 CN probes

n23 <- as.integer(sapply(cesN, FUN=getAttribute, "n23"));
cols <- c("blue", "red")[n23];
names(cols) <- colnames(theta);
cols["NA12145"] <- "purple";

units <- isSnp;
Mlim <- c(-2,2);
plotDensity(M[units,], col=cols, xlim=Mlim);


