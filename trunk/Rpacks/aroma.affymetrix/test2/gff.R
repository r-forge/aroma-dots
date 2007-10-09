if (interactive())
  savehistory();
library(aroma.affymetrix);
source("readGff.R");
source("appendCdfUnits.R");


key <- list(minLength=100e3, minCount=30);
dirs <- "CNVs";
tracks2 <- loadCache(key=key, dirs=dirs);

if (is.null(tracks2)) {
  if (!exists("cdf")) {
    cdf <- AffymetrixCdfFile$fromChipType("Mapping250K_Nsp");
  }
  
  pathname <- "500KEA_CNVs.gff";
  
  pathname <- Arguments$getReadablePathname(pathname, mustExist=TRUE);
  
  # Read all CNVs
  if (!exists("tracks")) {
    tracks <- readGff(pathname);
    rm(tracks);
  }
  
  # Extract units well within CNV regions
  if (!exists("tracks2")) {
    tracks2 <- appendCdfUnits(tracks, cdf, margin=100e3, minCount=30);
  }

  saveCache(tracks2, key=key, dirs=dirs);
}

nbrOfUnits <- sum(sapply(tracks2[[1]]$cnvs, .subset2, "nbrOfUnits"));
print(nbrOfUnits);

