if (interactive()) {
  savehistory();
} else {
  # GLAD v1.12.0 depends on the 'tcltk' package, but that
  # cannot be loaded if there is no display.  A workaround
  # is to fake that the 'tcltk' package is loaded:
  attach(list(), name="package:tcltk");
}
library(aroma.affymetrix);
#source("init.R");

# Use special file cache for testing
options("R.cache::rootPath"="~/.Rcache,scratch");
options("R.cache::touchOnLoad"=TRUE);


args <- commandArgs(asValues=TRUE, excludeReserved=TRUE, exludeEnvVars=TRUE);
print(args);

paths <- c();
allPaths <- c("testScripts/replication/chipTypes", "testScripts/system/chipTypes");
for (path in allPaths) {
  path <- Arguments$getReadablePath(path, mustExist=TRUE);
  paths0 <- list.files(path=path, full.names=TRUE);
  paths <- c(paths, paths0);
}

..pathnames <- lapply(paths, FUN=list.files, pattern="[.]R$", full.names=TRUE);
names(..pathnames) <- basename(paths);
..pathnames <- ..pathnames[names(..pathnames)];

..chipTypes <- c("Mapping10K_Xba142", "Test3",
                 "HG-U133_Plus_2", "Mapping50K_Hind240,Xba240",
                 "Mapping250K_Nsp,Sty", "HuEx-1_0-st-v2",
                 "GenomeWideSNP_6", "GenomeWideSNP_5");

..chipTypes <- rev(..chipTypes);

if (!is.null(args$chipTypes)) {
  ..chipTypes <- trim(unlist(strsplit(args$chipTypes, split=",")));
}

cat("Processing chip types:\n");
print(..chipTypes);

for (..chipType in ..chipTypes) {
  nbrOfTests <- length(..pathnames[[..chipType]]);
  for (kk in seq(length=nbrOfTests)) {
    pathname <- ..pathnames[[..chipType]][[kk]];
    if (regexpr("hetero", pathname) != -1)
      next;

    tryCatch({
      source(pathname, echo=TRUE);
    }, error = function(ex) {
      cat("************************************************\n");
      cat("** ", rep(c(" ER", "ROR "), times=6), " **\n", sep="");
      print(ex);
      cat("************************************************\n");
    });
    
    rm(list=setdiff(ls(), 
        c("kk", "..chipType", "..chipTypes", "pathname", "..pathnames")));
    gc();
  }
}
