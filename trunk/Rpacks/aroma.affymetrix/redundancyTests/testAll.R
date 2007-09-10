savehistory();
library(aroma.affymetrix);
source("init.R");

path <- "testScripts/system/chipTypes";
path <- Arguments$getReadablePath(path, mustExist=TRUE);
paths <- list.files(path=path, full.names=TRUE);

..pathnames <- lapply(paths, FUN=list.files, pattern="[.]R$", full.names=TRUE);
names(..pathnames) <- basename(paths);
..pathnames <- ..pathnames[names(..pathnames)];

print(names(..pathnames));

..chipTypes <- c("Test3", "Mapping10K_Xba142", "Mapping50K_Hind240,Xba240", "HG-U133_Plus_2", "GenomeWideSNP6.0");

pathname <- ..pathnames[[1]][1];
#source(pathname, echo=TRUE);
#stop();

..chipTypes <- ..chipTypes[3];
..chipTypes <- ..chipTypes[5];
..chipTypes <- ..chipTypes[1];

..chipTypes <- c("GenomeWideSNP6.0");
#..chipTypes <- c("Mapping10K_Xba142");

for (..chipType in ..chipTypes) {
  nbrOfTests <- length(..pathnames[[..chipType]]);
  for (kk in seq(length=nbrOfTests)) {
#    source("../aroma.affymetrix/R/digest2.R");
#    source("../aroma.affymetrix/R/AffymetrixCelSet.R");
    pathname <- ..pathnames[[..chipType]][[kk]];
    source(pathname, echo=TRUE);
    rm(list=setdiff(ls(), 
        c("kk", "..chipType", "..chipTypes", "pathname", "..pathnames")));
    gc();
  }
}
