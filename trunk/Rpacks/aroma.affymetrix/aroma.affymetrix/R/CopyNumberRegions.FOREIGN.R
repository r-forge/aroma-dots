setMethodS3("extractCopyNumberRegions", "profileCGH", function(object, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pv <- object$profileValues;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate result table
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify unique regions
  uRegions <- unique(pv$Region);
  nbrOfRegions <- length(uRegions);

  # Columns
  colClasses <- c(chromosome="character", start="integer", 
                  stop="integer", mean="double", nbrOfLoci="integer",
                  call="character");
  df <- dataFrame(colClasses, nrow=nbrOfRegions);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract each region
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (rr in seq(along=uRegions)) {
    # Get the region ID
    region <- uRegions[rr];

    # Get the first and last position of each region
    idx <- which(region == pv$Region);
    idx <- idx[c(1,length(idx))];
    idx1 <- idx[1];

    # Chromosome
    df[rr,"chromosome"] <- pv$Chromosome[idx1];

    # (start, stop, length)
    df[rr,c("start", "stop")] <- as.integer(pv$PosBase[idx]);

    # Number of SNPs
    df[rr,"nbrOfLoci"] <- as.integer(diff(idx)+1);

    # Smoothing
    df[rr,"mean"] <- pv$Smoothing[idx1];

    # Call
    df[rr,"call"] <- c("loss", "neutral", "gain")[pv$ZoneGNL[idx1]+2];
  }

  CopyNumberRegions(
    chromosome=df$chromosome,
    start=df$start, 
    stop=df$stop, 
    mean=df$mean, 
    count=df$nbrOfLoci,
    call=df$call
  );
})


setMethodS3("extractCopyNumberRegions", "DNAcopy", function(object, ...) {
  output <- object$output;

  CopyNumberRegions(
    chromosome=output[["chrom"]], 
    start=output[["loc.start"]], 
    stop=output[["loc.end"]], 
    mean=output[["seg.mean"]],
    count=output[["num.mark"]]
  );
})




############################################################################
# HISTORY:
# 2008-05-17
# o Extracted CopyNumberRegions.FOREIGN.R after moving CopyNumberRegions.R
#   to aroma.core.
# 2008-04-17
# o BUG FIX: applyRows() and lines() of CopyNumberRegions did not handle
#   cases with zero regions.
# 2007-09-04
# o Now CopyNumberRegions also contains an 'chromosome' field.
# o BUG FIX: as.data.frame() gave an error if some optional fields were
#   NULL.
# 2007-08-22
# o Created.  Need a generic container for holding copy number regions and
#   to plot them nicely.
############################################################################
