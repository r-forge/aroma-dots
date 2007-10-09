setMethodS3("appendCdfUnits", "default", function(tracks, cdf, margin=0, minCount=1, ...) {
  if (!inherits(cdf, "AffymetrixCdfFile")) {
    throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  gi <- getGenomeInformation(cdf);
  gp <- getData(gi);
  units <- as.integer(rownames(gp));

  # Add number of Nsp loci
  tracks <- lapply(tracks, FUN=function(track) { 
    data <- track$data;
    data[,"seqname"] <- gsub("^chr", "", data[,"seqname"]);
    data[,"seqname"] <- gsub("X", "23", data[,"seqname"]);
    data[,"seqname"] <- gsub("Y", "24", data[,"seqname"]);
    data[,"seqname"] <- as.integer(data[,"seqname"]);

    # For each CNV
    cnvs <- vector("list", nrow(data));

    # Process chromosome by chromosome 
    chromosomes <- sort(unique(data[,"seqname"]));
    for (chr in chromosomes) {
      rr <- which(data[,"seqname"] == chr);
      starts <- data[rr,"start"] + margin;
      ends <- data[rr,"end"] - margin;

      gpRows <- which(gp[,"chromosome"] == chr);
      pos <- gp[gpRows,"physicalPosition"];

      for (kk in seq(along=rr)) {
        keep <- which(starts[kk] <= pos & pos <= ends[kk]);
        cnvUnits <- units[gpRows[keep]];
        if (length(cnvUnits) >= minCount) {
          cnv <- list(
            chromosome=chr, 
            start=starts[kk],
            end=ends[kk], 
            length=ends[kk]-starts[kk]+1, 
            units=cnvUnits,
            nbrOfUnits=length(cnvUnits)
          );
        } else {
          cnv <- NA;
        }
        cnvs[[rr[kk]]] <- cnv;
      }
    } # for (chr in ...)

    cnvs <- cnvs[!is.na(cnvs)];

    track$cnvs <- cnvs;

    track;
  })

  tracks;
}) # appendCdfUnits()


############################################################################
# HISTORY:
# 2007-10-08
# o Created.  Will only be use for the private HapMap/CNV validation study.
############################################################################
