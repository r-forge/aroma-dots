library(aroma.affymetrix)
log <- Arguments$getVerbose(-4);
timestampOn(log);
.Machine$float.eps <- sqrt(.Machine$double.eps);

dataSetName <- "Jeremy_2007-10k";
chipType <- "Mapping10K_Xba142";

# Expected sample names
sampleNames <- c("0001-7", "0002-10", "0004-13", "0005-14", "0007-18", 
                      "0008-19", "0010-22", "2-DPrrr", "MH12", "MH18");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tests for setting up CEL sets and locating the CDF file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cs <- AffymetrixCelSet$byName(dataSetName, chipType=chipType, verbose=log);
keep <- 1:6;
cs <- extract(cs, keep);
sampleNames <- sampleNames[keep];
print(cs);
stopifnot(identical(getNames(cs), sampleNames));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spatial intensity plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ae <- ArrayExplorer(cs);
setColorMaps(ae, "sqrt,yellow");
print(ae);
stopifnot(identical(unname(getArrays(ae)), getNames(cs)));
process(ae, verbose=log);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spatial probe log-ratio plots
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pngDev <- findPngDevice();
cfR <- getAverageFile(cs, verbose=log);
reporter <- SpatialReporter(cs, reference=cfR);
ylab <- expression(log[2](y/y[R]));
figPath <- "figures";
for (array in 1:nbrOfArrays(cs)) {
  df <- getFile(cs, array);

  filename <- sprintf("%s,spatial,rowMedians.png", getFullName(df));
  pathname <- filePath(figPath, filename);
  devNew("pngDev", pathname, width=300, height=800);
  plotMargins(reporter, array=array, ylim=c(-1,1)*0.2, ylab=ylab, margins="rows", rotate=90);
  devDone();

  filename <- sprintf("%s,spatial,colMedians.png", getFullName(df));
  pathname <- filePath(figPath, filename);
  devNew("pngDev", pathname, width=800, height=300);
  plotMargins(reporter, array=array, ylim=c(-1,1)*0.2, ylab=ylab, margins="columns", rotate=0);
  devDone();
}

addColorMap(reporter, "log2center,rainbow");
process(reporter, zrange=c(-2,2), verbose=log);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Spatial residual plots test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plm <- RmaPlm(cs);
print(plm);
fit(plm, verbose=log);
rs <- calculateResidualSet(plm, verbose=log);
ae <- ArrayExplorer(rs);
setColorMaps(ae, c("log2,log2neg,rainbow", "log2,log2pos,rainbow"));
print(ae);
stopifnot(identical(unname(getArrays(ae)), getNames(cs)));
process(ae, interleaved="auto", verbose=log);
