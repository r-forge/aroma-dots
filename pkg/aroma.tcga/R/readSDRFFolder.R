readSDRFFolder <- function(allFiles=NULL, foldername, pattern="[_.,](sdrf|SDRF)[.](txt|TXT)$", verbose=FALSE) {
  if(is.null(allFiles)) {
    allFiles <- list.files(path=foldername, pattern=pattern, full.names=TRUE);
  }

  if(verbose) 
    cat("Reading Files...\n");

  tabList <- lapply(allFiles, FUN=readSDRFFile, verbose=verbose);

  tabCompile <- do.call(rbind, lapply(tabList, FUN=function(x){x}));

  if(verbose) 
    cat("Reading Files...done\n");

  tabCompile;
} # readSDRFFolder()

############################################################################
# HISTORY:
# 2008-03-19
# o Added argument 'pattern'.
# o Created from code by Elizabeth Purdom [EP].
############################################################################
