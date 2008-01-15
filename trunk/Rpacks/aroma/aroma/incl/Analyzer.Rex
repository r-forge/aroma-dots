# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# main()
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

if (!exists("ma")) {
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw, bg.subtract=TRUE)
  setTreatments(ma, c(1,1,1,2,2,2));
}

if (!exists("ss")) {
  ss  <- SpotSlideArray(ma$M, layout);
}

if (!exists("grs")) {
  grs <- as.GeneReplicateSlideArray(ss);
}

if (!exists("gs")) {
  gs <- as.GeneSlideArray(grs);
#  mat <- matrix(NA, nrow=length(grs), ncol=1);
#  tmp <- unclass(grs);
#  for (k in seq(length(grs))) {
#    mat[k,1] <- tmp[[k]][1];
#  }
#  gs <- GeneSlideArray(mat, layout);
}

if (!exists("an")) {
  an <- Analyzer(ma);
}

if (!exists("tma")) {
  tma <- ttest(an);
}


