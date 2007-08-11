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
setMethodS3("plotProfile2", "profileCGH", function(fit, variable="LogRatio", chromosome=NULL, Smoothing="Smoothing", GNL="ZoneGNL", Bkp=FALSE, cytobandLabels=TRUE, plotband=TRUE, unit=0, colDAGLAD=NULL, pchSymbol=c(20, 4), colCytoBand=c("white", "darkblue"), colCentro="red", xlim=NULL, ylim=c(-1,1)*2.5, xlab="Physical position", ylab=variable, flavor=c("glad", "ce", "minimal"), xmargin=c(50,50), resScale=1, ...) {
  require("GLAD") || stop("Package not loaded: GLAD");  # data(cytoband)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fit':
  if (!"PosBase" %in% names(fit$profileValues))
    throw("Argument 'fit' does not contain a 'PosBase' field.");

  # Argument 'variable':
  if (!variable %in% names(fit$profileValues))
    throw("Argument 'variable' does not specify a known field: ", variable);
  
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

  # Argument 'Smoothing':
  if (!is.null(Smoothing)) {
    if (!Smoothing %in% names(fit$profileValues)) {
      cat("Warning in plotProfile.profileCGH:", Smoothing, " is not available");
    }
  }

  # Argument 'Bkp':
  if (Bkp) {
    if (!"Breakpoints" %in% names(fit$profileValues))
      throw("Cannot plot breakpoints: No data available.")
  }

  # Argument 'colDAGLAD':
  if (is.null(colDAGLAD)) {
    colDAGLAD <- RColorBrewer::brewer.pal(5, "Dark2");
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Reset graphical parameters when done
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  opar <- par();
  on.exit(par(opar));


  xScale <- 1/(10^unit);

  # Keep only data to be plotted
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  pv <- fit$profileValues;

  # Keep only data for the chromosome of interest
  keep <- (pv$Chromosome == chromosome);
  pv <- pv[keep,];

  # Keep only finite values based on the variable of interest
  keep <- is.finite(pv[[variable]]);
  pv <- pv[keep,];

  # Convert the chromosome names to chromosome indices
  pv$Chromosome <- GLAD::ChrNumeric(pv$Chromosome);

  # Make sure the order of the values are increasing
  o <- order(pv$PosOrder);
  pv <- pv[o,];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get chromosome lengths
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Load data
  cytoband <- NULL; rm(cytoband); # To please R CMD check R v2.6.0 dev.
  data(cytoband);  # Package 'GLAD'
  genomeInfo <- aggregate(cytoband$End, list(Chromosome=cytoband$Chromosome, 
                          ChrNumeric=cytoband$ChrNumeric), max, na.rm=TRUE);
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
  Chromosome <- NULL; # dummy
  cytobandNew <- subset(cytoband, select=-Chromosome); 
##  cytobandNew <- cytoband[,setdiff(names(cytoband), "Chromosome"),drop=FALSE];

  cytobandNew <- merge(LabelChr, cytobandNew, by.x="Chromosome", 
                                                        by.y="ChrNumeric");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Update the plot data with cytoband information
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  pv <- merge(pv, LabelChr, by="Chromosome");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot smoothing values?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!is.null(Smoothing)) {
    # Order the data points by chromosome and position
    o <- order(pv$Chromosome, pv$PosBase);
    pv <- pv[o,];

    n <- length(pv[, 1])
    PosMax <- max(pv$PosBase) + 1;
    Pos <- pv$PosBase[1:(n-1)];
    PosNext <- pv$PosBase[2:n];
    InterPos <- Pos + (PosNext - Pos)/2;
    InterPos <- c(0, InterPos, PosMax);
    SmtStart <- pv[, Smoothing][1];
    SmtEnd <- pv[, Smoothing][n];
    Smt1 <- pv[, Smoothing][1:(n-1)];
    Smt1 <- c(SmtStart, Smt1, SmtEnd);
    Smt2 <- pv[, Smoothing][2:n];
    Smt2 <- c(SmtStart, Smt2, SmtEnd);
    datasmt <- data.frame(
      PosBase=c(InterPos, InterPos), 
      Smoothing=c(Smt1, Smt2)
    );
    rm(PosNext, SmtStart, SmtEnd, InterPos, Smt1, Smt2); # Not needed anymore
    datasmt <- unique(datasmt);
    datasmt <- datasmt[order(datasmt$PosBase),];
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Annotate gains, normals, and losses, as well as outliers?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (GNL %in% names(pv)) {
    # Setup color vector
    col <- rep(colDAGLAD[5], length(pv$PosOrder));
    gnl <- pv[GNL];
    col[gnl ==  -1] <- colDAGLAD[4];
    col[gnl ==   1] <- colDAGLAD[3];
    col[gnl ==   2] <- colDAGLAD[2];
    col[gnl == -10] <- colDAGLAD[1];

    # Setup pch vector
    pch <- rep(pchSymbol[1], length(pv$PosOrder));
    pch[pv$OutliersTot != 0] <- pchSymbol[2];

    colSmoothing <- "black";
  } else {
    col <- par("col");
    pch <- 20;
    colSmoothing <- "red";
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plotting flavor
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  opar <- par(no.readonly=TRUE);  
  on.exit(opar);

  if (flavor == "glad") {
    par(mar=c(3,3,5,3)+0.1, mgp=c(2,0.6,0.3));
    axes <- TRUE;
    if (is.null(xlim))
      xlim <- c(0, xScale*genomeInfo$Length[chromosome]);
  } else if (flavor == "ce") {
    # Margins in pixels-to-inches

    par(mar=c(3,3,5,3)+0.1, mgp=c(2,0.6,0.3), xaxs="i");

    # Set the horizontal margins to 'xmargin'.
    dim <- getDeviceResolution(resScale) * par("din");
    plt <- par("plt");    
    plt[1:2] <- c(xmargin[1], dim[1]-xmargin[2]) / dim[1];
    par("plt"=plt);

    axes <- TRUE;
    if (is.null(xlim))
      xlim <- c(0, xScale*genomeInfo$Length[chromosome]);
  } else if (flavor == "minimal") {
    # No margins
    	par(mar=c(2,0,0.5,0), mgp=c(2,0.6,0.3), xaxs="i");
    # No axis
    axes <- FALSE;
    # No cytobands
    plotband <- FALSE;
    # x-range
    if (is.null(xlim))
      xlim <- c(0, xScale*genomeInfo$Length[chromosome]);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot main data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the data to plot
  y <- pv[, variable];
  x <- xScale*pv$PosBase;
  plot(x=x, y=y, pch=pch, col=col, xaxt="n", xlab=xlab, ylab=ylab, 
#                            xlim=xlim, ylim=ylim, axes=axes, bty="n", ...);
                            xlim=xlim, ylim=ylim, axes=axes, ...);
#  axis(side=2); axis(side=4);
   

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot cytobands?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot break points?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (Bkp) {
    if (is.data.frame(fit$BkpInfo)) {
      fit$BkpInfo <- merge(fit$BkpInfo, LabelChr, by="Chromosome");
      fit$BkpInfo$NewPosBase <- fit$BkpInfo$PosBase + fit$BkpInfo$Length;
      fit$BkpInfo$NewPosBase <- xScale*fit$BkpInfo$NewPosBase;
      abline(v=fit$BkpInfo$NewPosBase + 0.5, col="red", lty=2);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot smoothing lines?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!is.null(Smoothing)) {
    datasmt$PosBase <- xScale * datasmt$PosBase;
    lines(datasmt$Smoothing ~ datasmt$PosBase, col=colSmoothing);
  }
}, private=TRUE) # plotProfile2()


############################################################################
# HISTORY:
# 2007-06-11
# o Added explicit call to GLAD::myPalette() to please R CMD check R v2.6.0.
# 2007-01-03
# o Made the highlighting "arrow" for the centromere smaller.
# 2006-12-20
# o It is now possible to specify 'xlim' as well as 'ylim'.
# o Reimplemented, because the cytoband was not displayed correctly.
############################################################################
