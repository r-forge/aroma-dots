pruneSDRF <- function(sdrfAll, samplefile, todo=c("keep", "delete"), colSdrf="Array.Data.File", colSampleFile=colSdrf) {
  todo <- match.arg(todo);

  sdrfMatch <- which(sdrfAll[,colSdrf] %in% samplefile[,colSampleFile]);
  if(length(sdrfMatch) == 0) 
    stop("No matches found");

  if(todo == "keep") 
    sdrfKeep <- sdrfAll[sdrfMatch,];

  if(todo=="delete") 
    sdrfKeep <- sdrfAll[-sdrfMatch,];

  sdrfKeep;
} # pruneSDRF()


############################################################################
# HISTORY:
# 2008-03-19
# o Created from code by Elizabeth Purdom [EP].
############################################################################
