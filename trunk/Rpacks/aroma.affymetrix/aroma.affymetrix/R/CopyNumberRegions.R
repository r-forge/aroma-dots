setConstructorS3("CopyNumberRegions", function(start=NULL, stop=NULL, mean=NULL, count=NULL, call=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  extend(Object(), "CopyNumberRegions", 
    start = start,
    stop = stop,
    mean = mean,
    count = count,
    call= call,
    ...
  )
})


setMethodS3("as.character", "CopyNumberRegions", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Number of regions: %d", nbrOfRegions(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE) 


setMethodS3("nbrOfRegions", "CopyNumberRegions", function(this, ...) {
  length(this$start);
})


setMethodS3("as.data.frame", "CopyNumberRegions", function(x, ...) {
  # To please R CMD check
  this <- x;

  data <- cbind(start=this$start, stop=this$stop, mean=this$mean, count=this$count, call=this$call);
  data;
})


setMethodS3("applyRegions", "CopyNumberRegions", function(this, FUN, ...) {
  data <- as.data.frame(this);
  o <- order(data[,"start"]);
  data <- data[o,];
  apply(data, MARGIN=1, FUN=FUN);
})


setMethodS3("drawLevels", "CopyNumberRegions", function(this, col="red", lwd=2, lty=1, ...) {
  col0 <- col;
  lwd0 <- lwd;
  lty0 <- lty;
  applyRegions(this, FUN=function(cnr) {
    x <- c(cnr["start"], cnr["stop"]);
    y <- rep(cnr["mean"], times=2);
    if (is.function(col0))
      col <- col0(cnr);
    if (is.function(lwd0))
      lwd <- lwd0(cnr);
    if (is.function(lty0))
      lty <- lty0(cnr);
    lines(x=x, y=y, col=col, lwd=lwd, lty=lty, ...);
  })
})



setMethodS3("lines", "CopyNumberRegions", function(x, col="red", lwd=2, ...) {
  # To please R CMD check.
  this <- x;

  data <- as.data.frame(this);
  o <- order(data[,"start"]);
  data <- data[o,];
  xx <- t(data[,c("start", "stop")]);
  yy <- rep(this$mean[o], each=2);
  lines(x=xx, y=yy, col=col, lwd=lwd, ...);
})



setMethodS3("extractCNRs", "default", function(...) {
  extractCopyNumberRegions(...);
})


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
    start=output[["loc.start"]], 
    stop=output[["loc.end"]], 
    mean=output[["seg.mean"]],
    count=output[["num.mark"]]
  );
})




############################################################################
# HISTORY:
# 2007-08-22
# o Created.  Need a generic container for holding copy number regions and
#   to plot them nicely.
############################################################################
