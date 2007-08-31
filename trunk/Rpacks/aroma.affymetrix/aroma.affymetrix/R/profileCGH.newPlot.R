setMethodS3("drawCytoBand", "default", function(cytoband, chromosome=1, y=-1, labels=TRUE, height=1, colCytoBand=c("white", "darkblue"), colCentro="red", ...) {
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


# Patch for plotProfile() of class profileCGH so that 'ylim' argument works.
# Added also par(cex=0.8) - see code.
setMethodS3("newPlot", "profileCGH", function(fit, chromosome=NULL, GNL="ZoneGNL", cytobandLabels=TRUE, plotband=TRUE, unit=0, colCytoBand=c("white", "darkblue"), colCentro="red", xlim=NULL, ylim=c(-1,1)*2.5, xlab="Physical position", ylab="Relative copy number", flavor=c("glad", "ce", "minimal"), xmargin=c(50,50), resScale=1, ...) {
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

  # Argument 'flavor':
  flavor <- match.arg(flavor);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Reset graphical parameters when done
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  opar <- par();
  on.exit(par(opar));


  xScale <- 1/(10^unit);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get chromosome lengths
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Load data
  # To please R CMD check on R v2.6.0
  cytoband <- NULL; rm(cytoband);
  data("cytoband");  # Package 'GLAD'
  genomeInfo <- aggregate(cytoband$End, list(Chromosome=cytoband$Chromosome, 
                          ChrNumeric=cytoband$ChrNumeric), max, na.rm=TRUE);
  names(genomeInfo) <- c("Chromosome", "ChrNumeric", "Length");
  genomeInfo$Chromosome <- as.character(genomeInfo$Chromosome);
  genomeInfo$ChrNumeric <- as.integer(as.character(genomeInfo$ChrNumeric));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plotting flavor
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  opar <- par(no.readonly=TRUE);  
  on.exit(opar);

  # x-range
  if (is.null(xlim)) {
    xlim <- c(0, xScale*genomeInfo$Length[chromosome]);
  }

  if (flavor == "glad") {
    par(mar=c(3,3,5,3)+0.1, mgp=c(2,0.6,0.3));
    axes <- TRUE;
  } else if (flavor == "ce") {
    # Margins in pixels-to-inches

    par(mar=c(3,3,5,3)+0.1, mgp=c(2,0.6,0.3), xaxs="i");

    # Set the horizontal margins to 'xmargin'.
    dim <- getDeviceResolution(resScale) * par("din");
    plt <- par("plt");    
    plt[1:2] <- c(xmargin[1], dim[1]-xmargin[2]) / dim[1];
    par("plt"=plt);

    axes <- TRUE;
  } else if (flavor == "minimal") {
    # No margins
    	par(mar=c(2,0,0.5,0), mgp=c(2,0.6,0.3), xaxs="i");
    # No axis
    axes <- FALSE;
    # No cytobands
    plotband <- FALSE;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot main data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the data to plot
  plot.window(xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                                                 axes=axes, xaxt="n", ...);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the cytoband details for the chromosome of interest
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Drop column 'Chromosome'
  ## Gives a NOTE in R CMD check R v2.6.0, which is nothing, but we'll
  ## use a workaround to get a clean result. /HB 2007-06-12
  Chromosome <- NULL; rm(Chromosome); # dummy
  cytobandNew <- subset(cytoband, select=-Chromosome); 
##  cytobandNew <- cytoband[,setdiff(names(cytoband), "Chromosome"),drop=FALSE];

  LabelChr <- data.frame(Chromosome=chromosome);
  LabelChr <- merge(LabelChr, genomeInfo[, c("ChrNumeric", "Length")], 
                         by.x="Chromosome", by.y="ChrNumeric", all.x=TRUE);

  LabelChr$Length <- 0;


  cytobandNew <- merge(LabelChr, cytobandNew, by.x="Chromosome", 
                                                        by.y="ChrNumeric");
  if (plotband) {
    # Rescale x positions according to units
    cytobandNew$Start <- xScale*cytobandNew$Start;
    cytobandNew$End <- xScale*cytobandNew$End;

    usr <- par("usr");
    dy <- diff(usr[3:4]);

    drawCytoBand(cytobandNew, chromosome=chromosome, 
      labels=cytobandLabels, y=usr[4]+0.02*dy, height=0.03*dy, 
      colCytoBand=colCytoBand, colCentro=colCentro);
  }
}, private=TRUE) # newPlot()


############################################################################
# HISTORY:
# 2007-08-23
# o Created from profileCGH.plotProfile2.R.
############################################################################
