
setMethodS3("calcMvsA", "AffymetrixCelFile", function(this, reference, indices=NULL, ..., zeros=FALSE) {
  # Arguments 'reference':
  if (!inherits(reference, "AffymetrixCelFile")) {
    throw("Argument 'reference' is not of class AffymetrixCelFile: ", 
                                                          class(reference)[1]);
  }

  # Argument 'indices':
  nbrOfCells <- nbrOfCells(this);
  if (is.null(indices)) {
  } else {
    indices <- Arguments$getIndices(indices, range=c(1,nbrOfCells));
  }

  # Check if the two CEL files are compatible
  if (nbrOfCells != nbrOfCells(reference)) {
    throw("This and the 'reference' CEL file have different number of cells: ", 
                                   nbrOfCells, " != ", nbrOfCells(reference));
  }

  # Get the signals for this channel
  y1 <- getData(this, indices=indices, fields="intensities")[,1];

  offset <- this$offset;
  if (is.null(offset))
    offset <- 0;
  if (offset != 0)
    cat("Offset: ", offset, "\n", sep="");

  # Identify non-zero signals?
  if (!zeros) {
    keep <- which(y1 != 0);
    y1 <- y1[keep];
  } else {
    keep <- seq(along=y1);
  }
  y1 <- y1 + offset;
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
    y2 <- getData(reference, indices=indices, fields="intensities")[,1];
    y2 <- y2 + offset;
    y2 <- log(y2, base=2);
  }

  ma <- matrix(c((y1+y2)/2, y1-y2), ncol=2);
  colnames(ma) <- c("A", "M");
  ma;
})

setMethodS3("annotateMvsA", "AffymetrixCelFile", function(this, reference, ..., what="M") {
  if (identical(what, "M")) {
    abline(h=0, lty=1, col="blue");
  }
  stextChipType(this);
  stextLabels(this, others=reference);
}, protected=TRUE)

setMethodS3("plotMvsA", "AffymetrixCelFile", function(this, reference, indices=NULL, pch=176, xlim=c(0,16), ylim=c(-1,1)*diff(xlim), xlab=expression(A==1/2%*%log[2](y[1]*y[2])), ylab=expression(M==log[2](y[1]/y[2])), ..., annotate=TRUE) {
  ma <- calcMvsA(this, reference, indices=indices);
  plot(ma, pch=pch, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...);
  if (annotate)
    annotateMvsA(this, reference);
  this$lastPlotData <- ma;
  invisible(ma);
})

setMethodS3("smoothScatterMvsA", "AffymetrixCelFile", function(this, reference, indices=NULL, pch=176, xlim=c(0,16), ylim=c(-1,1)*diff(xlim), xlab=expression(A==1/2%*%log[2](y[1]*y[2])), ylab=expression(M==log[2](y[1]/y[2])), ..., annotate=TRUE) {
  require(geneplotter) || throw("Package 'geneplotter' not loaded.");
  ma <- calcMvsA(this, reference, indices=indices);
  smoothScatter(ma, pch=pch, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...);
  if (annotate)
    annotateMvsA(this, reference);
  this$lastPlotData <- ma;
  invisible(ma);
})


setMethodS3("stextChipType", "AffymetrixCelFile", function(this, side=4, fmtstr="%s", pos=1, cex=0.7, col="darkgray", ...) {
  stext(side=side, text=sprintf(fmtstr, getChipType(this)), pos=pos, cex=cex, col=col, ...);
})


setMethodS3("plotMvsX", "AffymetrixCelFile", function(this, reference, x, indices=NULL, pch=176, ylim=c(-1,1)*2, ylab=expression(M==log[2](y1/y2)), ..., what="M", add=FALSE, annotate=!add) {
  ma <- calcMvsA(this, reference, indices=indices, zeros=TRUE);
  nobs <- nrow(ma);
  if (nobs == 0)
    throw("Cannot plot M vs X because there is not non-zero data.");

  if (add) {
    points(x, ma[,what], pch=pch, ...);
  } else {
    plot(x, ma[,what], pch=pch, ylim=ylim, ylab=ylab, ...);
    if (annotate)
      annotateMvsA(this, reference, what=what);
  }

  # The first two columns should always be the data plotted
  ma <- cbind(x=x, M=ma[,what], A=ma[,setdiff(colnames(ma),what)]);

  this$lastPlotData <- ma;
  invisible(ma);
})

setMethodS3("plotMvsPosition", "AffymetrixCelFile", function(this, reference, chromosome, region=NULL, gdas, ylab=expression(M==log[2](y1/y2)), xlab="Physical position [Mb]", xlim=NULL, ..., what="M", annotate=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'region':
  if (!is.null(region)) {
    if (is.character(region)) {
    } else {
      region <- Arguments$getDoubles(region, range=c(0,Inf), length=c(2,2));
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get GDAS annotations for the specific chromosome.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  df <- select(gdas, Chromosome=chromosome);

  # Get the loci (in Mb)
  loci <- as.integer(df[,"PhysicalPosition"]);
  loci <- loci / 10^6;

  # Record genome size
  chrRange <- range(loci, na.rm=TRUE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Filter out by region of interest?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(region)) {
    # Region specified by SNP names?
    if (is.character(region)) {
      # One, two, or several SNPs may be given
      # Find matching SNP names (ignore missing)
      rr <- match(region, rownames(df));
      rr <- na.omit(rr);
      if (length(rr) == 0) {
        throw("Argument 'region' contains unknown SNPs for chromosome '", chromosome, "': ", paste(region, collapse=", "));
      }
      region <- range(loci[rr]);

      # If only one SNP was specified/identified, show 10% of the chromosome
      if (diff(region) == 0)
        region <- region + c(-1,1)*0.05*chrRange;
    }

    keep <- (region[1] <= loci & loci <= region[2]);
    loci <- loci[keep];
    df <- df[keep,,drop=FALSE];

    if (is.null(xlim))
      xlim <- region;
  } else {
    if (is.null(xlim))
      xlim <- c(0,chrRange[2]);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Re-order by physical position
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  o <- order(loci);
  loci <- loci[o];
  unitNames <- rownames(df)[o];
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove non-existing data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Make sure all units exists in data too
  cdf <- getCdf(this);
  allUnitNames <- getUnitNames(cdf);
  units <- match(unitNames, allUnitNames);

  # Just in case there are no existing units
  keep <- which(!is.na(units));
  units <- units[keep];
  loci <- loci[keep];
  unitNames <- unitNames[keep];

  if (length(units) == 0)
    throw("Could not identify units for this request.");
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the cell indices for these
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cells <- unlist(getCellIndices(this, units=units), use.names=FALSE);
  if (length(cells) == 0)
    throw("Could not identify cell indices for this request.");

  # Assert correctness of data dimensions before continuing
  if (length(cells) != length(loci)) {
    throw("Internal error: Cannot plot M vs genomic position. The number of cells does not match the number of loci: ", length(cells), " != ", length(loci));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot data against location
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- plotMvsX(this, reference, x=loci, indices=cells, ylab=ylab, xlab=xlab, xlim=xlim, ..., what=what, annotate=annotate);
  if (annotate) {
    stext(side=3, pos=1, sprintf("Chromosome %s", chromosome), cex=0.7);
    x <- res[,what];
    mad <- mad(x[is.finite(x)]);
    stext(sprintf("MAD: %.3g", mad), side=3, line=-1, pos=1, 
                                                  cex=0.7, col="darkgray");
  }

  rownames(res) <- unitNames;
  this$lastPlotData <- res;

  invisible(res);
})


setMethodS3("highlight", "AffymetrixCelFile", function(this, indices=NULL, ...) {
  data <- this$lastPlotData;
  if (!is.null(indices))
    data <- data[indices,,drop=FALSE];
  points(data[,1:2], ...);
  invisible(data);
})

# match(c("SNP_A-1737080", "SNP_A-1740686"), rownames(df))
############################################################################
# HISTORY:
# 2006-08-27
# o Added plotMvsX() and plotMvsPosition().
# o Added calcMvsA(), plotMvsA(), and smoothScatterMvsA().
############################################################################
