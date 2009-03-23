# Details:
# How is the Regional GC Correction algorithm enhancement implemented?
# The Regional GC Correction enhancement in the copy number algorithm utilizes information on GC content from the NetAffx NA26.1 annotation. For each marker, the annotation file reports the percent GC content in a half-megabase window. The CN algorithm bins the markers based on this GC content and normalizes the log2 ratios within each bin. Then the algorithm calculates an adjustment factor and applies this factor to each log2 ratio. This generates the GC-corrected log2 ratio for each marker. Please see the Copy Number Algorithm GC Waviness Correction white paper for complete details on the GC Correction algorithm enhancement.
# Reference:
# http://www.affymetrix.com/support/help/faqs/genotyping_console/copy_number_analysis/faq_5.jsp

setMethodS3("plotCovariateEffects", "GcContentNormalization2", function(this, arrays=NULL, units=NULL, ref="zero", ..., xlim=NULL, ylim=NULL, xlab="GC fraction", ylab=expression(log[2](theta)), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dataSet <- getInputDataSet(this);
  unf <- getUnitNamesFile(dataSet);
  # Argument 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(unf)));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Plotting signals as a function of covariate(s)");

  X <- getCovariates(gcn, units=units, verbose=less(verbose,5));
  verbose & cat(verbose, "Covariates:");
  verbose & str(verbose, X);
  X <- X[,1];

  keep <- is.finite(X);
  X <- X[keep];
  units <- units[keep];
  rm(keep);

  yR <- NULL;
  if (ref == "median") {
    verbose && enter(verbose, "Extracting reference signals");
    ceR <- getAverageFile(dataSet, verbose=less(verbose,10));
    yR <- extractTotalAndFreqB(ceR, units=units, drop=TRUE, verbose=less(verbose,50))[,"total"];
    yR <- log2(yR);
    verbose && str(verbose, yR);
    verbose && exit(verbose);
    if (is.null(ylim)) {
      ylim <- c(-1,1)*3;
    }
  } else if (ref == "zero") {
    if (is.null(ylim)) {
      ylim <- c(0,16);
    }
  }

  if (is.null(xlim)) {
    xlim <- range(X, na.rm=TRUE);
  }

  verbose && enter(verbose, "Plotting");
  nbrOfArrays <- nbrOfArrays(dataSet);
  subplots(nbrOfArrays);
  for (cc in seq(length=nbrOfArrays)) {
    ce <- getFile(dataSet, cc);
    name <- getFullName(ce);
    name <- gsub(",chipEffects", "", name);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", cc, name, nbrOfArrays));
    x <- X;

    verbose && enter(verbose, "Extracting total signals");
    y <- extractTotalAndFreqB(ce, units=units, drop=TRUE, verbose=less(verbose,50))[,"total"];
    y <- log2(y);
    if (!is.null(yR)) {
      y <- y - yR;
    }
    verbose && exit(verbose);
    ok <- (is.finite(x) & is.finite(y));
    x <- x[ok];
    y <- y[ok];
    rm(ok);
    verbose && cat(verbose, "(x,y):");
    verbose && str(verbose, x);
    verbose && str(verbose, y);
    geneplotter::smoothScatter(x, y, pch=".", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim);
    stext(side=3, pos=1, name);
    fit <- smooth.spline(x,y);
    lines(fit, col="red", lwd=2);
    rm(x,y);
    verbose && exit(verbose);
  } # for (cc ...)
  verbose && exit(verbose);
  rm(yR);

  verbose && exit(verbose);
}) # plotCovariateEffects()


############################################################################
# HISTORY:
# 2009-03-22
# o Created.
############################################################################
