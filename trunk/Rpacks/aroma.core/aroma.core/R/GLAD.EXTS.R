setMethodS3("drawCnRegions", "profileCGH", function(this, xScale=1, ...) {
  # Get data
  pv <- this$profileValues;

  # Order data along chromosome(s)
  o <- order(pv$Chromosome, pv$PosBase);
  pv <- pv[o, ];

  # Number of data points
  n <- length(pv[, 1]);

  posMax <- max(pv$PosBase) + 1;
  pos <- pv$PosBase[1:(n-1)];
  posNext <- pv$PosBase[2:n];

  interPos <- pos + (posNext - pos)/2;
  interPos <- c(0, interPos, posMax);
  smtStart <- pv[, "Smoothing"][1];
  smtEnd <- pv[, "Smoothing"][n];
  smt1 <- pv[, "Smoothing"][1:(n-1)];
  smt1 <- c(smtStart, smt1, smtEnd);
  smt2 <- pv[, "Smoothing"][2:n];
  smt2 <- c(smtStart, smt2, smtEnd);
  datasmt <- data.frame(posBase=c(interPos, interPos), smoothing=c(smt1, smt2));
  datasmt <- unique(datasmt);
  datasmt <- datasmt[order(datasmt$posBase), ];
#  posBase <- xScale * datasmt$posBase;

  lines(datasmt$smoothing ~ datasmt$posBase, ...);
})



# Patch for plotProfile() of class profileCGH so that 'ylim' argument works.
# Added also par(cex=0.8) - see code.
setMethodS3("drawCytoband", "profileCGH", function(fit, chromosome=NULL, cytobandLabels=TRUE, colCytoBand=c("white", "darkblue"), colCentro="red", unit=6, ...) {
  require("GLAD") || stop("Package not loaded: GLAD");  # data("cytoband")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fit':
  if (!"PosBase" %in% names(fit$profileValues))
    throw("Argument 'fit' does not contain a 'PosBase' field.");

  # Argument 'chromosome':
  if (is.null(chromosome)) {
    chromosome <- unique(fit$profileValues$Chromosome);
    if (length(chromosome) > 1) {
      throw("Argument 'chromosome' must not be NULL if 'fit' contains more than one chromosome: ", paste(chromosome, collapse=", "));
    }
  }
  if (length(chromosome) > 1) {
    throw("Argument 'chromosome' must not contain more than one chromosome: ", paste(chromosome, collapse=", "));
  }


  xScale <- 1/(10^unit);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get chromosome lengths
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Load data
  # To please R CMD check on R v2.6.0
  cytoband <- NULL; rm(cytoband);
  data("cytoband", envir=sys.frame(sys.nframe()));  # Package 'GLAD'
  genomeInfo <- aggregate(cytoband$End, 
    by=list(Chromosome=cytoband$Chromosome, ChrNumeric=cytoband$ChrNumeric), 
    FUN=max, na.rm=TRUE);
  names(genomeInfo) <- c("Chromosome", "ChrNumeric", "Length");
  genomeInfo$Chromosome <- as.character(genomeInfo$Chromosome);
  genomeInfo$ChrNumeric <- as.integer(as.character(genomeInfo$ChrNumeric));

  LabelChr <- data.frame(Chromosome=chromosome);
  LabelChr <- merge(LabelChr, genomeInfo[, c("ChrNumeric", "Length")], 
                         by.x="Chromosome", by.y="ChrNumeric", all.x=TRUE);

  LabelChr$Length <- 0;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the cytoband details for the chromosome of interest
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Drop column 'Chromosome'
  ## Gives a NOTE in R CMD check R v2.6.0, which is nothing, but we'll
  ## use a workaround to get a clean result. /HB 2007-06-12
  Chromosome <- NULL; rm(Chromosome); # dummy
  cytobandNew <- subset(cytoband, select=-Chromosome); 
  cytobandNew <- merge(LabelChr, cytobandNew, by.x="Chromosome", 
                                                        by.y="ChrNumeric");
  # Rescale x positions according to units
  cytobandNew$Start <- xScale*cytobandNew$Start;
  cytobandNew$End <- xScale*cytobandNew$End;

  # Where should the cytoband be added and how wide should it be?
  usr <- par("usr");
  dy <- diff(usr[3:4]);

  drawCytoband2(cytobandNew, chromosome=chromosome, 
    labels=cytobandLabels, y=usr[4]+0.02*dy, height=0.03*dy, 
    colCytoBand=colCytoBand, colCentro=colCentro);
}, private=TRUE) # drawCytoband()



setMethodS3("drawCytoband2", "default", function(cytoband, chromosome=1, y=-1, labels=TRUE, height=1, colCytoBand=c("white", "darkblue"), colCentro="red", ...) {
  opar <- par(xpd=NA);
  on.exit(par(opar));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Cytoband colors
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  color <- unique(cytoband$Col);
  pal <- GLAD::myPalette(low=colCytoBand[1], high=colCytoBand[2], k=length(color));

  info <- data.frame(Color=color, ColorName=I(pal));
  cytoband <- merge(cytoband, info, by="Color");
  rm(info);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract cytoband information for current chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  keep <- which(cytoband$Chromosome == chromosome);
  cytoband <- cytoband[keep, ];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Cytoband positions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  CytoPos <- 0.5 * (cytoband$Start + cytoband$End);
  CytoLength <- (cytoband$End - cytoband$Start);
  NbCyto <- length(cytoband[, 1]);
  HeightPlot <- rep(height, NbCyto);
  sizeCyto <- matrix(c(CytoLength, HeightPlot), nrow=NbCyto, ncol=2);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Draw cytobands
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  y0 <- min(unique(y));
  yC <- y0+height/2;
  y1 <- y0+height;
  symbols(x=CytoPos, y=rep(yC, NbCyto), rectangles=sizeCyto,
      inches=FALSE, bg=cytoband$ColorName, add=TRUE, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Highlight the centromere
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # The inverted arrow indicating where the centromere is.
  idxs <- which(cytoband$Centro == 1);
  centroPos <- min(cytoband$End[idxs]);
  arrows(centroPos, y0, centroPos, y1, col=colCentro, code=2, angle=120, 
                                                               length=0.1);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Labels, e.g. 20q12
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (labels) {
    labels <- paste(cytoband$Chromosome, cytoband$Band, sep="");
#    axis(side=3, at=CytoPos, labels=labels, las=2);
    dy <- par("cxy")[2];
    text(x=CytoPos, y=y1+dy/2, labels=labels, srt=90, adj=c(0,0.5));
  }
}, private=TRUE)

 
############################################################################
# HISTORY:
# 2009-05-10
# o Moved to aroma.core v1.0.6.  Source files: profileCGH.drawCnRegions.R
#   and profileCGH.drawCytoband.R.
# 2007-09-04
# o Now data("cytoband") is loaded to the local environment.
# 2007-08-22
# o Update plotProfile2() to utilizes drawCnRegions().
# 2007-08-20
# o Added drawCnRegions().
# 2007-06-11
# o Added explicit call to GLAD::myPalette() to please R CMD check R v2.6.0.
# 2007-01-03
# o Made the highlighting "arrow" for the centromere smaller.
# 2006-12-20
# o It is now possible to specify 'xlim' as well as 'ylim'.
# o Reimplemented, because the cytoband was not displayed correctly.
############################################################################ 
