###########################################################################/**
# @set "class=AffymetrixCelFile"
# @RdocMethod plotDensity
#
# @title "Plots the density of the probe signals on the array"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{subset}{The subset of probes to include.}
#   \item{types}{The type of probes to include.}
#   \item{...}{Additional arguments passed to 
#      \code{identifyCells()} of @see "AffymetrixCdfFile",
#      , @seemethod "getData", and 
#      @see "aroma.light::plotDensity.numeric".}
#   \item{xlab,ylab}{The labels on the x and the y axes.}
#   \item{log}{If @TRUE, the density of the log (base 2) values are 
#      used, otherwise the non-logged values.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("plotDensity", "AffymetrixCelFile", function(this, subset=1/2, types=NULL, ..., xlab=NULL, ylab="density (integrates to one)", log=TRUE, annotate=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'subset':

  # Argument 'xlab':
  if (is.null(xlab)) {
    if (log) {
      xlab <- expression(log[2](y));
    } else {
      xlab <- expression(y);
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);
  verbose && enter(verbose, "Identifying subset of probes");
  suppressWarnings({
    subset <- identifyCells(cdf, indices=subset, types=types, ...,
                                                    verbose=less(verbose));
  })
  verbose && exit(verbose);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot density
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Plotting the density");
  verbose && cat(verbose, "Array: ", getName(this));
  suppressWarnings({
    verbose && enter(verbose, "Loading probe intensities");
    y <- getData(this, indices=subset, fields="intensities");
    y <- y$intensities;
    verbose && exit(verbose);
    if (log) {
      verbose && cat(verbose, "Taking the logarithm (base 2)");
      y <- log(y, base=2);
    }
    verbose && cat(verbose, "Plotting");
    plotDensity(y, xlab=xlab, ylab=ylab, ...);
  })

  if (annotate) {
    stextChipType(getCdf(this));
    stextLabels(this);
    stextSize(this, size=length(y));
  }

  verbose && exit(verbose);
})



setMethodS3("getAm", "AffymetrixCelFile", function(this, reference, indices=NULL, ..., zeros=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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

#  # Argument 'indices' & 'units':
#  if (!is.null(indices) && !is.null(units)) {
#    throw("Arguments 'indices' and 'units' must not be non-NULL at the same time.");
#  }
#
#  # Argument 'units':
#  cdf <- getCdf(this);
#  nbrOfUnits <- nbrOfUnits(cdf);
#  if (is.null(units)) {
#  } else {
#    units <- Arguments$getIndices(units, range=c(1,nbrOfUnits));
#  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Further validation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Check if the two CEL files are compatible
  if (nbrOfCells != nbrOfCells(reference)) {
    throw("This and the 'reference' CEL file have different number of cells: ", 
                                   nbrOfCells, " != ", nbrOfCells(reference));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the signals for this channel
  y1 <- getData(this, indices=indices, fields="intensities")[,1];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Offset signals?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  offset <- this$offset;
  if (is.null(offset))
    offset <- 0;
  if (offset != 0)
    cat("Offset: ", offset, "\n", sep="");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Remove signals that are zero?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!zeros) {
    keep <- which(y1 != 0);
    y1 <- y1[keep];
  } else {
    keep <- seq(along=y1);
  }
  y1 <- y1 + offset;
  y1 <- log(y1, base=2);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get reference signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return (A,M)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  am <- matrix(c((y1+y2)/2, y1-y2), ncol=2);
  colnames(am) <- c("A", "M");

  am;
})



setMethodS3("annotateMvsA", "AffymetrixCelFile", function(this, reference, ..., what="M") {
  if (identical(what, "M")) {
    abline(h=0, lty=1, col="blue");
  }
  stextChipType(getCdf(this));
  stextLabels(this, others=reference);
}, private=TRUE)




###########################################################################/**
# @RdocMethod plotMvsA
#
# @title "Plots log-ratio versus log-intensity in a scatter plot"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{reference}{The reference channel, i.e. the denominator in the
#     log ratios.}
#   \item{indices}{Indices of the probe signals to be plotted.}
#   \item{pch}{The plot symbol.}
#   \item{xlim,ylim}{The range of the x and the y axes.}
#   \item{xlab,ylab}{The labels on the x and the y axes.}
#   \item{...}{Additional arguments passed to @see "graphics::plot".} 
#   \item{annotate}{If @TRUE, the plot is annotated with information about
#     the data plotted, otherwise not.}
# }
#
# \value{
#  Returns (invisibly) a @data.frame with the plotted data in the
#  first two columns.
# }
#
# @author
#
# \seealso{
#   @seemethod "smoothScatterMvsA".
#   @seemethod "plotMvsX".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("plotMvsA", "AffymetrixCelFile", function(this, reference, indices=NULL, pch=176, xlim=c(0,16), ylim=c(-1,1)*diff(xlim), xlab=expression(A==1/2%*%log[2](y[1]*y[2])), ylab=expression(M==log[2](y[1]/y[2])), ..., annotate=TRUE) {
  ma <- getAm(this, reference, indices=indices);
  plot(ma, pch=pch, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...);
  if (annotate) {
    annotateMvsA(this, reference);
    stextSize(this, size=nrow(ma));
  }
  this$lastPlotData <- ma;
  invisible(ma);
})



###########################################################################/**
# @RdocMethod smoothScatterMvsA
#
# @title "Plots log-ratio versus log-intensity in a smooth scatter plot"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{reference}{The reference channel, i.e. the denominator in the
#     log ratios.}
#   \item{indices}{Indices of the probe signals to be plotted.}
#   \item{pch}{The plot symbol.}
#   \item{xlim,ylim}{The range of the x and the y axes.}
#   \item{xlab,ylab}{The labels on the x and the y axes.}
#   \item{...}{Additional arguments passed to @see "graphics::plot".} 
#   \item{annotate}{If @TRUE, the plot is annotated with information about
#     the data plotted, otherwise not.}
# }
#
# \value{
#  Returns (invisibly) a @data.frame with the plotted data in the
#  first two columns.
# }
#
# @author
#
# \seealso{
#   @seemethod "plotMvsA".
#   Internally @see "geneplotter::smoothScatter" is used.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("smoothScatterMvsA", "AffymetrixCelFile", function(this, reference, indices=NULL, pch=176, xlim=c(0,16), ylim=c(-1,1)*diff(xlim), xlab=expression(A==1/2%*%log[2](y[1]*y[2])), ylab=expression(M==log[2](y[1]/y[2])), ..., annotate=TRUE) {
  require(geneplotter) || throw("Package 'geneplotter' not loaded.");
  ma <- getAm(this, reference, indices=indices);
  smoothScatter(ma, pch=pch, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...);
  if (annotate) {
    annotateMvsA(this, reference);
    stextSize(this, size=nrow(ma));
  }
  this$lastPlotData <- ma;
  invisible(ma);
})




###########################################################################/**
# @RdocMethod plotMvsX
#
# @title "Plots log-ratio versus another variable in a scatter plot"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{reference}{The reference channel, i.e. the denominator in the
#     log ratios.}
#   \item{x}{The other variable.  A @double @vector.}
#   \item{indices}{Indices of the probe signals to be plotted.}
#   \item{pch}{The plot symbol.}
#   \item{xlim,ylim}{The range of the x and the y axes.}
#   \item{xlab,ylab}{The labels on the x and the y axes.}
#   \item{...}{Additional arguments passed to @see "graphics::plot".} 
#   \item{annotate}{If @TRUE, the plot is annotated with information about
#     the data plotted, otherwise not.}
# }
#
# \value{
#  Returns (invisibly) a @data.frame with the plotted data in the
#  first two columns, and remaining data in the following columns.
# }
#
# @author
#
# \seealso{
#   @seemethod "plotMvsA".
#   @seemethod "smoothScatterMvsA".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("plotMvsX", "AffymetrixCelFile", function(this, reference, x, indices=NULL, pch=176, ylim=c(-1,1)*2, ylab=NULL, ..., what=c("M", "A"), add=FALSE, annotate=!add) {
  # Argument 'what':
  what <- match.arg(what);

  # Get the log-ratios
  ma <- getAm(this, reference, indices=indices, zeros=TRUE);
  nobs <- nrow(ma);
  if (nobs == 0)
    throw("Cannot plot M vs X because there is not non-zero data.");

  if (nobs != length(x)) {
    throw("The number of log-ratios does not match the number of elements in argument 'x': ", length(nobs), " != ", length(x));
  }

  if (what == "M") {
    ylab <- expression(M==log[2](y1/y2))
  } else {
    ma <- ma[,2:1];
    ylab <- expression(A==1/2%*%log[2](y1*y2))
  }

  if (add) {
    points(x, ma[,1], pch=pch, ...);
  } else {
    plot(x, ma[,1], pch=pch, ylim=ylim, ylab=ylab, ...);
    if (annotate) {
      annotateMvsA(this, reference, what=what);
      stextSize(this, size=length(x));
    }
  }

  # The first two columns should always be the data plotted
  ma <- cbind(x=x, ma);

  this$lastPlotData <- ma;
  invisible(ma);
})



###########################################################################/**
# @RdocMethod plotMvsPosition
#
# @title "Plots log-ratio versus physical position in a scatter plot"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{reference}{The reference channel, i.e. the denominator in the
#     log ratios.}
#   \item{x}{The other variable.  A @double @vector.}
#   \item{chromosome}{The chromosome to be plotted.}
#   \item{region}{The region to be plotted.  If @NULL, the whole chromosome
#     is plotted.}
#   \item{pch}{The plot symbol.}
#   \item{xlim,ylim}{The range of the x and the y axes.}
#   \item{xlab,ylab}{The labels on the x and the y axes.}
#   \item{...}{Additional arguments passed to @see "graphics::plot".} 
#   \item{annotate}{If @TRUE, the plot is annotated with information about
#     the data plotted, otherwise not.}
# }
#
# \value{
#  Returns (invisibly) a @data.frame with the plotted data in the
#  first two columns, and remaining data in the following columns.
# }
#
# @author
#
# \seealso{
#   @seemethod "plotMvsX".
#   @seemethod "plotMvsA".
#   @seemethod "smoothScatterMvsA".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("plotMvsPosition", "AffymetrixCelFile", function(this, reference, chromosome, region=NULL, gdas, ylab=ylab, xlab="Physical position [Mb]", xlim=NULL, ..., what="M", annotate=TRUE) {
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

  # Argument 'what':
  what <- match.arg(what);

  if (what == "M") {
    ylab <- expression(M==log[2](y1/y2))
  } else {
    ma <- ma[,2:1];
    ylab <- expression(A==1/2%*%log[2](y1*y2))
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
    stextSize(this, size=length(cells));
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



############################################################################
# HISTORY:
# 2006-09-26
# o Renamed calcMvsA() to getAm().
# 2006-09-15
# o Added more Rdoc comments.
# o Readded plotDensity(). 
# o Added stextSize() to annotate with "n=1034".
# 2006-08-27
# o Added plotMvsX() and plotMvsPosition().
# o Added calcMvsA(), plotMvsA(), and smoothScatterMvsA().
# 2006-07-27
# o Added argument 'verbose' to plotDensity().
# 2006-05-29
# o Added Rdoc comments.
# 2006-05-16
# o Created.
############################################################################
