readSDRFFile <- function(filename, verbose=FALSE) {
  if(verbose)
    cat(filename,"\n");

  tab <- read.table(filename, stringsAsFactors=FALSE, header=TRUE, na.strings="->", sep="\t", fill=TRUE);

  if("Protocol.REF.2" %in% names(tab)) {
    tab <- tab[, c("Extract.Name", "Hybridization.Name", "Array.Data.File", "Protocol.REF.2")];
    names(tab)[4] <- "Protocol.REF"; # Instead? /HB
  } else {
    tab <- tab[, c("Extract.Name", "Hybridization.Name", "Array.Data.File", "Protocol.REF")];
  }

  # remove if don't start with TCGA
  tab <- tab[grep("^TCGA", tab$Extract.Name), ];


  df <- sapply(tab$Extract.Name, FUN=function(name) {
    unlist(parseSdrfExtractName(name));
  });
  colnames(df) <- NULL;
  df <- t(df);

##  center <- sapply(strsplit(tab$Extract.Name, split=""), FUN=function(x) {
##    paste(tail(x, n=2), collapse="");
##  });
##
##  # first take off the center stuff
##  Portion.id <- sub("[A-Z]-[0-9]*-[0-9]{2}$", "", tab$Extract.Name);
##  Sample.id <- sub("-[0-9]{2}$", "", Portion.id);
##  Patient.id <- sub("-[0-9]{2}[A-Z]$", "", Sample.id);

  center.Name <- sapply(strsplit(tab$Protocol.REF, ":"), FUN=.subset2, 1);
  type <- sapply(strsplit(tab$Protocol.REF, ":"), FUN=.subset2, 3);
  type <- sub("_[0-9]$", "", type);

  tab <- data.frame(tab, 
    Center.id=df[,"centerId"], Center.name=center.Name, Patient.id=df[,"patientId"], 
    Sample.id=df[,"sampleId"], Portion.id=df[,"portionId"], Platform.Type=type, 
    stringsAsFactors=FALSE
  );

  tab <- tab[, c("Patient.id", "Hybridization.Name", "Platform.Type", 
                 "Protocol.REF", "Array.Data.File", "Extract.Name", 
                 "Center.id", "Center.name", "Sample.id", "Portion.id")];

  tab;
} # readSDRFFile()


############################################################################
# HISTORY:
# 2008-03-19
# o BUG FIX: Added fill=TRUE to read.table() [HB].
# o Created from code by Elizabeth Purdom [EP].
############################################################################
