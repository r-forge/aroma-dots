setMethodS3("calcMvsA", "AffymetrixCelFile", function(this, reference, indices=NULL, ...) {
  # Arguments 'reference':
  if (!inherits(reference, "AffymetrixCelFile")) {
    throw("Argument 'reference' is not of class AffymetrixCelFile: ", 
                                                          class(reference)[1]);
  }

  # Check if the two CEL files are compatible
  if (nbrOfCells(this) != nbrOfCells(reference)) {
    throw("This and the 'reference' CEL file have different number of cells: ", 
                              nbrOfCells(this), " != ", nbrOfCells(reference));
  }

  # Get the signals for this channel
  y1 <- getFields(this, indices=indices, fields="intensities")[,1];

  # Identify non-zero signals
  keep <- which(y1 != 0);
  y1 <- y1[keep];
  y1 <- log(y1, base=2);

  if (length(y1) == 0) {
    y2 <- y1;
  } else {
    # Get the signals for the reference channel
    if (is.null(indices)) {
      indices <- keep;
    } else {
      indices <- indices[keep];
    }
    y2 <- getFields(reference, indices=indices, fields="intensities")[,1];
    y2 <- log(y2, base=2);
  }

  ma <- matrix(c((y1+y2)/2, y1-y2), ncol=2);
  colnames(ma) <- c("A", "M");
  ma;
})

setMethodS3("annotateMvsA", "AffymetrixCelFile", function(this, reference, ...) {
  abline(h=0, lty=1, col="blue");
  stextChipType(this);
  stextLabels(this, others=reference);
}, protected=TRUE)

setMethodS3("plotMvsA", "AffymetrixCelFile", function(this, reference, indices=NULL, pch=176, xlim=c(0,16), ylim=c(-1,1)*diff(xlim), xlab=expression(A==1/2*log[2](y[1]*y[2])), ylab=expression(M==log[2](y[1]/y[2])), ..., annotate=TRUE) {
  ma <- calcMvsA(this, reference, indices=indices);
  plot(ma, pch=pch, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...);
  if (annotate)
    annotateMvsA(this, reference);
})

setMethodS3("smoothScatterMvsA", "AffymetrixCelFile", function(this, reference, indices=NULL, pch=176, xlim=c(0,16), ylim=c(-1,1)*diff(xlim), xlab=expression(A==1/2*log[2](y[1]*y[2])), ylab=expression(M==log[2](y[1]/y[2])), ..., annotate=TRUE) {
  require(geneplotter) || throw("Package 'geneplotter' not loaded.");
  ma <- calcMvsA(this, reference, indices=indices);
  smoothScatter(ma, pch=pch, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...);
  if (annotate)
    annotateMvsA(this, reference);
})


setMethodS3("stextChipType", "AffymetrixCelFile", function(this, side=4, fmtstr="Chip type: %s", pos=1, cex=0.7, col="darkgray", ...) {
  stext(side=side, text=sprintf(fmtstr, getChipType(this)), pos=pos, cex=cex, col=col, ...);
})


setMethodS3("plotMvsX", "AffymetrixCelFile", function(this, reference, x, indices=NULL, ylim=c(-1,1)*2, ylab=expression(M==log[2](y1/y2)), ..., annotate=TRUE) {
  ma <- calcMvsA(this, reference, indices=indices);
  plot(x, ma[,"M"], ylim=ylim, ylab=ylab, ...);
  if (annotate)
    annotateMvsA(this, reference);
})



setMethodS3("plotMvsPosition", "AffymetrixCelFile", function(this, reference, chromosome, gdas, ylab=expression(M==log[2](y1/y2)), xlab="Physical position [Mb]", ...) {
  dfAll <- gdas[];

  rr <- which(dfAll[,"Chromosome"] == as.character(chromosome));
  df <- dfAll[rr,];

  # Get the loci
  loci <- as.integer(df[,"PhysicalPosition"]);
  o <- order(loci);
  loci <- loci[o];

  unitNames <- rownames(df)[o];
  
  # Make sure all units exists in data too
  cdf <- getCdf(this);
  allUnitNames <- getUnitNames(cdf);
  units <- match(unitNames, allUnitNames);
  keep <- which(!is.na(units));
  units <- units[keep];
  loci <- loci[keep];
  loci <- loci / 10^6; # In Mb

  # Identify the cell indices for these
  cells <- unlist(getFirstCellIndices(ce, units=units), use.names=FALSE);

  plotMvsX(this, reference, x=loci, indices=cells, ylab=ylab, xlab=xlab, ...);
})


############################################################################
# HISTORY:
# 2006-08-27
# o Added plotMvsX() and plotMvsPosition().
# o Added calcMvsA(), plotMvsA(), and smoothScatterMvsA().
############################################################################
