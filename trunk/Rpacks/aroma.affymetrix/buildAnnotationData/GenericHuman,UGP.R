#############################################################################
# Author: Henrik Bengtsson, hb@stat.berkeley.edu
# Created: 2008-05-20
# Description: This script generated UGP files that are generic to any 
# platform and that contains unit locations that are uniformly spread
# along the genome at various resolutions, e.g. on locus every 50kb.
#############################################################################
if (interactive()) savehistory();
library("aroma.core");

log <- verbose <- Arguments$getVerbose(-20, timestamp=TRUE);

chipType <- "GenericHuman";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Table containing the lengths of each chromosomes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("annotationData", "genomes", "Human");
path <- Arguments$getReadablePath(path);  
filename <- "Human,chromosomes,max,20090503.txt";
hg <- TabularTextFile(filename, path=path);
df <- readDataFrame(hg);
chr <- as.character(df[,"chromosome"]);
chr[chr == "X"] <- 23;
chr[chr == "Y"] <- 24;
chr[chr == "M"] <- 25;
df[,"chromosome"] <- as.integer(chr);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("annotationData", "chipTypes", chipType);
path <- Arguments$getWritablePath(path);  

# One unit every X kb in [1,#bases] starting at base 1.
for (by in c(5,10,50,100,200)*1e3) {
  ugpData <- NULL;
  for (kk in seq(along=df[,"chromosome"])) {
    chr <- df[kk,"chromosome"];
    pos <- seq(from=1, to=df[kk,"nbrOfBases"], by=by);
    ugpKK <- data.frame(chromosome=chr, pos=pos);
    ugpData <- rbind(ugpData, ugpKK);
  }

  chipType <- sprintf("GenericHuman,%dkb", as.integer(by/1e3)); 
  filename <- sprintf("%s,HB20090503.ugp", chipType);
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
  if (!isFile(pathname)) {
    footer <- list(createdOn=format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE), createdBy="Henrik Bengtsson, hb@stat.berkeley.edu", genome="Human", chromosomes=paste(sort(unique(df[,"chromosome"])), collapse=","), description=sprintf("A (unit, chromosome, position) table where units are located every %dkb across the whole human genome", as.integer(by/1e3)));
    ugp <- AromaUgpFile$allocate(filename=pathname, path=NULL, nbrOfRows=nrow(ugpData), platform="Generic", chipType=chipType, footer=footer, overwrite=TRUE, verbose=log);

    ugp[,1] <- ugpData[,1];
    ugp[,2] <- ugpData[,2];
    print(ugp);
  }
} # for (by ...)

#############################################################################
# HISTORY:
# 2009-05-03
# o Added loci on Chr25 (mitochondrial DNA).
# 2008-05-20
# o Created.
#############################################################################
