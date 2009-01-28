if (interactive()) savehistory();
library("aroma.cn");

platforms <- list(
  "A" = list(
    dataSet  = "broad.mit.edu_GBM.Genome_Wide_SNP_6.1.1.0,summarized",
    chipType = "GenomeWideSNP_6"
  ),

  "B" = list(
     dataSet = "mskcc.org_GBM.HG-CGH-244A.1.data",
    chipType = "HG-CGH-244A"
  ),

  "C" = list(
     dataSet = "stanford.edu_GBM.HumanHap550.2.Paired_LogR",
    chipType = "HumanHap550"
  )
)


for (pp in seq(along=platforms)) {
  platform <- platforms[[pp]];
  dataSet <- platform$dataSet;
  chipType <- platform$chipType;
  ds <- AromaUnitTotalCnBinarySet$byName(dataSet, chipType=chipType);
  print(ds);
  platform$ds <- ds;
  platforms[[pp]] <- platform;
}

#freqBSet <- AromaUnitFracBCnBinarySet$byName(dataSet, chipType=chipType);
#print(freqBSet);


ylim <- NULL;
layout(matrix(1:length(platforms), ncol=1));
for (pp in seq(along=platforms)) {
  platform <- platforms[[pp]];
  ds <- platform$ds;
  ugp <- getAromaUgpFile(ds);
  print(ugp);
  
  units <- getUnitsAt(ugp, 7);
  pos <- ugp[units,2, drop=TRUE];
  y <- extractMatrix(ds, units=units, files=1);
  plot(pos, y, pch=".", cex=2, ylim=ylim);
}
