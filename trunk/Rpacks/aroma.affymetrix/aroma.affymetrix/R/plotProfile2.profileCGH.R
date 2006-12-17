require(GLAD) || stop("Package GLAD not found");

# Patch for plotProfile() of class profileCGH so that 'ylim' argument works.
# Added also par(cex=0.8) - see code.
setMethodS3("plotProfile2", "profileCGH", function(fit, variable="LogRatio", chromosome=NULL, Smoothing=NULL, GNL="ZoneGNL", Bkp=FALSE, labels=TRUE, plotband=TRUE, unit=0, colDAGLAD=c("black", "blue", "red", "green", "yellow"), pchSymbol=c(20, 13), colCytoBand=c("white", "darkblue"), colCentro="red", text=NULL, main="", ylim=NULL, ...) {
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Keep only data to be plotted
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Keep only data for the chromosome of interest
  keep <- (fit$profileValues$Chromosome == chromosome);
  fit$profileValues <- fit$profileValues[keep,];

  # Extract the data to plot
  fit$profileValues$VarToPlot <- fit$profileValues[, variable];

  # Convert the chromosome names to chromosome indices
  fit$profileValues$Chromosome <- ChrNumeric(fit$profileValues$Chromosome);

  # Keep only finite values
  keep <- is.finite(fit$profileValues[variable]);
  fit$profileValues <- fit$profileValues[keep,];

  # Make sure the order of the values are increasing
  o <- order(fit$profileValues$PosOrder);
  fit$profileValues <- fit$profileValues[o,];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Cytoband data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Load data
  data(cytoband);

  LabelChr <- data.frame(Chromosome=chromosome);
  genomeInfo <- aggregate(cytoband$End, list(Chromosome=cytoband$Chromosome, ChrNumeric=cytoband$ChrNumeric), max, na.rm=TRUE);
  names(genomeInfo) <- c("Chromosome", "ChrNumeric", "Length");
  genomeInfo$Chromosome <- as.character(genomeInfo$Chromosome);
  genomeInfo$ChrNumeric <- as.integer(as.character(genomeInfo$ChrNumeric));

  LabelChr <- merge(LabelChr, genomeInfo[, c("ChrNumeric", "Length")], by.x="Chromosome", by.y="ChrNumeric", all.x=TRUE);

  LabelChr$Length <- 0;

  cytobandNew <- subset(cytoband, select=-Chromosome);
  cytobandNew <- merge(LabelChr, cytobandNew, by.x="Chromosome", by.y="ChrNumeric");

  fit$profileValues <- merge(fit$profileValues, LabelChr, by="Chromosome");
  fit$profileValues$NewPosBase <- fit$profileValues$PosBase + fit$profileValues$Length;
  par(no.readonly=TRUE);

  # Plot cytobands
  if (plotband) {
    # Rescale x positions according to units
    cytobandNew$Start <- cytobandNew$Start/(10^unit);
    cytobandNew$End <- cytobandNew$End/(10^unit);
    cytobandNew$Start <- cytobandNew$Start + cytobandNew$Length;
    cytobandNew$End <- cytobandNew$End + cytobandNew$Length;
  
    # Plot cytoband and main graph in two panels
    layout(c(1,2), heights=c(1,4));
    par(mar=c(0,4,4,2));
    plot(0, type="n", xlim=c(0, max(cytobandNew$End)), ylim=c(-1.5, 1.5), xaxt="n", yaxt="n", ylab="", xlab="");
    LabelChrCyto <- unique(cytobandNew$Chromosome);
opar <- par(cex=0.8); # HB
    plotCytoBand(cytobandNew, Chromosome=LabelChrCyto[1], 
      labels=labels, y=0, height=2, colCytoBand=colCytoBand, 
      colCentro=colCentro);
par(opar); #HB 
    par(mar=c(4,4,0,2));
  }

  if (!is.null(Smoothing)) {
    fit$profileValues <- fit$profileValues[order(fit$profileValues$Chromosome, fit$profileValues$PosBase),];
    NbPos <- length(fit$profileValues[, 1])
    PosMax <- max(fit$profileValues$NewPosBase) + 1;
    Pos <- fit$profileValues$NewPosBase[1:(NbPos-1)];
    PosNext <- fit$profileValues$NewPosBase[2:NbPos];
    InterPos <- Pos + (PosNext - Pos)/2;
    InterPos <- c(0, InterPos, PosMax);
    SmtStart <- fit$profileValues[, Smoothing][1];
    SmtEnd <- fit$profileValues[, Smoothing][NbPos];
    Smt1 <- fit$profileValues[, Smoothing][1:(NbPos-1)];
    Smt1 <- c(SmtStart, Smt1, SmtEnd);
    Smt2 <- fit$profileValues[, Smoothing][2:NbPos];
    Smt2 <- c(SmtStart, Smt2, SmtEnd);
    datasmt <- data.frame(PosBase=c(InterPos, InterPos), Smoothing=c(Smt1, Smt2));
    datasmt <- unique(datasmt);
    datasmt <- datasmt[order(datasmt$PosBase),];
  }

  if (GNL %in% names(fit$profileValues)) {
    col <- rep(colDAGLAD[5], length(fit$profileValues$PosOrder));
    gnl <- fit$profileValues[GNL];
    col[gnl ==  -1] <- colDAGLAD[4];
    col[gnl ==   1] <- colDAGLAD[3];
    col[gnl ==   2] <- colDAGLAD[2];
    col[gnl == -10] <- colDAGLAD[1];
    outliers <- rep(pchSymbol[1], length(fit$profileValues$PosOrder));
    outliers[fit$profileValues$OutliersTot != 0] <- pchSymbol[2];

    if (plotband) {
      plot(VarToPlot ~ NewPosBase, data=fit$profileValues, pch=outliers, col=col, xaxt="n", xlab=main, ylab=variable, ylim=ylim, ...);
    } else {
      plot(VarToPlot ~ NewPosBase, data=fit$profileValues, pch=outliers, col=col, xaxt="n", xlab="", ylab=variable, ylim=ylim, main=main);
    }

    if (!is.null(Smoothing)) {
      lines(datasmt$Smoothing ~ datasmt$PosBase, col="black");
    }

    if (Bkp) {
      if (is.data.frame(fit$BkpInfo)) {
        fit$BkpInfo <- merge(fit$BkpInfo, LabelChr, by="Chromosome");
        fit$BkpInfo$NewPosBase <- fit$BkpInfo$PosBase + fit$BkpInfo$Length;
        abline(v=fit$BkpInfo$NewPosBase + 0.5, col="red", lty=2);
      }
    }
  } else {
    if (plotband) {
      plot(VarToPlot ~ NewPosBase, data=fit$profileValues, pch=20, xaxt="n", xlab=main, ylab=variable, ylim=ylim, ...);
    } else {
      plot(VarToPlot ~ NewPosBase, data=fit$profileValues, pch=20, xaxt="n", xlab="", ylab=variable, ylim=ylim, main=main, ...);
    }

    if (Bkp) {
      if (is.data.frame(fit$BkpInfo)) {
        fit$BkpInfo <- merge(fit$BkpInfo, LabelChr, by="Chromosome");
        fit$BkpInfo$NewPosBase <- fit$BkpInfo$PosBase + fit$BkpInfo$Length;
        abline(v=fit$BkpInfo$NewPosBase + 0.5, col="red", lty=2);
      }
    }
    if (!is.null(Smoothing)) {
      lines(datasmt$Smoothing ~ datasmt$PosBase, col="red");
    }
  }
  if (!is.null(text)) {
    text(text$x, text$y, labels=text$labels, cex=text$cex);
  }
}) # plotProfile2()

