############################################################################
# This file contains:
#
#  1) normalizePlatewise()
#  2) normalizeGenewise()
#  3) normalizePrintorder() (almost like 1)
#  4) normalizeSpatial()
#
#  Still to come:
#  5) normalizePrinttipwise()
#
#  Currently method number 1 is making use of the help functions:
#
#     i) normalizeGroupsConstant()
#    ii) normalizeGroupsShiftFunction()
#
#  Later also methods 2, 3 & 5 will make use of these.
############################################################################

############################################################################
# normalizeGroupsConstant()
#
# This is a private general method for normalizing the data 'y' under
# the assumption that it has a constant shift and constant scale within
# each group. This function will both adjust the bias (if not FALSE) and
# the scale (if not FALSE).
############################################################################
setMethodS3("normalizeGroupsConstant", "MicroarrayData", function(this, y, spots, bias=TRUE, scale=TRUE, robust=TRUE) {
  ngroups <- length(spots);
  # Remove the bias group by group, and at the same time
  # calculate the discrepancy for each group...
  b <- c();
  s <- c();
  for (k in 1:ngroups) {
    idx <- spots[[k]];
    yg <- y[idx];
    ok <- !is.na(yg);
    if (robust == TRUE)
      c <- median(yg[ok])
    else
      c <- mean(yg[ok]);

    # Calculate the discrepancy for the current group.
    if (robust == TRUE)
      d <- 1.4826*median(abs(yg[ok]-c))
    else
      d <- var(yg[ok]-c);
    s <- c(s, d);
    
    if (bias == FALSE) {
      b <- c(b, c); # Still bias in data
    } else if (bias == TRUE || bias == 0) {
      # Remove the bias within the plate. Assumes zero shift is the truth.
      yg <- yg - c;
      y[idx] <- yg;
      b <- c(b, 0); # No bias left
    } else {
      # Remove the bias within the plate. Assumes zero shift is the truth.
      yg <- yg - c;
      # Add the wanted bias
      yg <- yg + bias;
      y[idx] <- yg;
      b <- c(b, bias); # No bias left
    }
  } # for (k in 1:ngroups)

  if (scale != FALSE) {
    # The geometrical mean of the plate discrepancies.
    s.wanted <- prod(s)^(1/length(s));

    # If a special scale is wanted, adjust...
    if (is.numeric(scale))
      s.wanted <- scale / s.wanted;
      
    # Rescale each plate to get scale equal to s.gmean...
    for (k in 1:ngroups) {
      idx <- spots[[k]];
      y[idx] <- b[k] + (y[idx] - b[k]) * (s.wanted / s[k]);
    }
  } # if (scale...)

  # Return the normalized data
  y;
}, private=TRUE, static=TRUE)   # normalizeGroupConstant()






############################################################################
# normalizeGroupsShiftFunction()
#
# This is a private general method for normalizing the data 'y' under the
# assumption that it has a shift which is a function of 'x', e.g. the
# intensity. This method will adjust the shift to have zero shift and
# all groups having the same scale.
############################################################################
setMethodS3("normalizeGroupsShiftFunction", "MicroarrayData", function(this, x, y, spots, span=0.3, scale=NULL, robust=TRUE) {
  ngroups <- length(spots);
  # Remove the bias group by group, and at the same time
  # calculate the discrepancy for each group...
  s <- c();
  for (k in 1:ngroups) {
    idx <- spots[[k]];
    xg <- x[idx];
    yg <- y[idx];
    ok <- !is.na(xg) & !is.infinite(xg) & !is.na(yg) & !is.infinite(yg);

    # Estimate the function y = c(x) + error
    c <- lowess(x=xg[ok], y=yg[ok], f=span);

    # Find the value of c(x) for each (complete) data point (x,y)
    yg.c <- approx(c, xout=xg[ok], ties=mean)$y;

    # Adjust the bias
    yg[ok] <- yg[ok]-yg.c;

    if (is.null(scale) || !is.na(scale)) {
      # Calculate the scale for the current group which is assumed to be constant.
      if (robust == TRUE)
        d <- 1.4826*median(abs(yg[ok]))
      else
        d <- var(yg[ok]);
      s <- c(s, d);
    }
    
    y[idx] <- yg;
  } # for (k in 1:ngroups)

  # Now when we know the scale of each group, give them the same scale.
  if (is.null(scale) || !is.na(scale)) {
    # The geometrical mean of the plate discrepancies.
    s.geom <- prod(s)^(1/length(s));

    if (is.numeric(scale))
      s.geom <- scale / s.geom;

    # Rescale each plate to get scale equal to s.gmean...
    for (k in 1:ngroups) {
      idx <- spots[[k]];
      y[idx] <- y[idx] * (s.geom / s[k]);
    }
  } # if (scale...)
  
  # Return the normalized data
  y;
}, private=TRUE, static=TRUE)   # normalizeGroupShiftFunction()




############################################################################
# normalizeGenewise()
############################################################################
setMethodS3("normalizeGenewise", "MicroarrayData", function(this, fields=getFieldNames(this), bias=0, scale=1, robust=TRUE, verbose=TRUE, ...) {
  # Assert that argument 'field' is correctly specified.
  if (length(fields) == 0) {
    warning("No field(s) to normalize.");
    return(NULL);
  } else {
    match <- match(fields, getFieldNames(this));
    nok <- is.na(match);
    if (any(nok))
      throw("Argument 'fields' contains unknown fields: ", fields[nok])
  }

  # Assert that argument 'bias' is correctly specified.
  if (!is.numeric(bias) & !is.na(bias))
    throw("Argument 'bias' must be numeric: ", bias)
  
  # Assert that argument 'scale' is correctly specified.
  if (!is.numeric(scale) & !is.na(scale))
    throw("Argument 'scale' must be numeric: ", scale)

  # Set the bias/scale for all fields.
  biases <- rep(bias, length.out=length(fields))
  scales <- rep(scale, length.out=length(fields))

  ###########################################################
  #
  ###########################################################
  layout <- getLayout(this);
  geneGroups <- getGeneGroups(layout);
  genes <- getSpots(geneGroups);
  ngenes <- nbrOfGroups(geneGroups);
  
  # Normalize field by field...
  for (k in seq(fields)) {
    field <- fields[k];
    if (verbose == TRUE)
      cat("Normalize field ", field, " genewise...", sep="");
    data  <- this[[field]];
    bias  <- biases[k];
    scale <- scales[k];
    # ...and gene by gene...
    if (robust == TRUE) {
      # Robust standardization
      for (l in 1:ngenes) {
        spots <- genes[[l]];
        x <- data[spots,];
        ok <- !is.na(x);
        center <- median(x[ok])
        x <- x - center;
        scale <- 1.4826 * median(abs(x[ok]));
        data[spots,] <- x/scale;
      }
    } else {
      # Non-robust standardization
      for (l in 1:ngenes) {
        spots <- genes[[l]];
        x <- data[spots,];
        ok <- !is.na(x);
        center <- median(x[ok])
        x <- x - center;
        scale <- 1.4826 * median(abs(x[ok]));
        data[spots,] <- x/scale;
      }
    }
    this[[field]] <- data;
    if (verbose == TRUE)
      cat("done\n");
  }

  clearCache(this); 

  return(invisible(this));
}, private=TRUE, trial=TRUE) # normalizeGenewise()



############################################################################
# normalizePrintorder()
############################################################################
# jumppoints is a vector of jump points, which specifies
# the segements to be normalized (independently).
setMethodS3("normalizePrintorder", "MicroarrayData", function(this, fields, breakpoints, slides=NULL, method="loess", bias=TRUE, scale=TRUE, meanValue=NULL, smooth=0.30, robust=TRUE, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  
  if (method == "loess") {
    require(modreg); # loess()
    if (robust == TRUE)
      family <- "symmetric"
    else
      family <- "gaussian";
  } else if (method == "ssanova") {
    require(modreg); # smooth.spline()
    throw("Spatial normalization using splines is not impleted: ", method);
  } else {
    throw("Unknown method for spatial normalization: ", method);
  }

  if (length(fields) == 0)
    throw("Missing field(s) to normalize.");

  if (is.null(meanValue)) {
  } else if (is.numeric(meanValue)) {
    meanValue <- rep(meanValue, length.out=length(fields));
  } else {
    throw("Invalid value of argument 'meanValue': ", meanValue);
  }
  meanValue0 <- meanValue;
  
  if (length(breakpoints) == 0)
    throw("Missing breakpoints.");

  breakpoints <- c(1, breakpoints, 1+nbrOfSpots(this));
  breakpoints <- unique(breakpoints);
  breakpoints <- sort(breakpoints);
  nsections <- length(breakpoints)-1;

  layout <- getLayout(this);
  printorder <- toPrintorderMatrix(layout);

  res <- list();
  if (bias == TRUE) {
    var <- NULL;
    for (slide in slides) {
      fieldIdx <- 0;
      for (field in fields) {
        fieldIdx <- fieldIdx + 1;
        data <- this[[field]][printorder,slide];
        fitted <- rep(NA, nbrOfSpots(this));
        for (k in 1:nsections) {
          idx <- breakpoints[k]:(breakpoints[k+1]-1);
          y <- data[idx];
          if (method == "loess") {
            ok <- !is.na(y);
            idx <- idx[ok];
            y <- data[idx];
            fit <- loess(y ~ idx, family=family, span=smooth);
            yres <- as.vector(residuals(fit));
            yfit <- fitted(fit);
            data[idx] <- yres;
            fitted[printorder[idx]] <- yfit;
            rm(yres); rm(yfit);
          }
        } # for (k in 1:nsections)

        # Adjust for mean values
        meanValue <- meanValue0;
        if (is.null(meanValue[fieldIdx]))
          meanValue[fieldIdx] <- mean(fitted, na.rm=TRUE);
        str(meanValue[fieldIdx]);
        data <- data + meanValue[fieldIdx];

        # Store the results
        this[[field]][printorder,slide] <- data;
        res[[field]] <- cbind(res[[field]], fitted);
        dimnames(res[[field]]) <- NULL;
      } # for (field in fields)
    } # for (slide in slides)
  }
  
  if (scale == TRUE) {
    var <- NULL;
    for (slide in slides) {
      for (field in fields) {
        data <- this[[field]][printorder,slide];
        for (k in 1:nsections) {
          idx <- breakpoints[k]:(breakpoints[k+1]-1);
          y <- data[idx];
          var <- c(var, mad(y, na.rm=TRUE));
        }
        Z <- prod(var)^(1/nsections);  # (a scalar)
        scale <- Z / var;
        for (k in 1:nsections) {
          idx <- breakpoints[k]:(breakpoints[k+1]-1);
          data[idx] <- scale[k] * data[idx];
        }
      } # for (field in fields)
    } # for (slide in slides)
  }

  clearCache(this); 

  invisible(res);
})
 

############################################################################
# normalizeSpatial()
############################################################################
setMethodS3("normalizeSpatial", "MicroarrayData", function(this, fields, slides=NULL, bias=NULL, method="loess", smooth=0.02, ...) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (method == "loess") {
    require(modreg); # loess()
  } else if (method == "ssanova") {
    require(gss); # ssanova()
    throw("Spatial normalization using thin plate splines (ssanova) is not impleted: ", method);
  } else {
    throw("Unknown value on argument 'method': ", method);
  }

  if (length(fields) == 0)
    throw("Missing field(s) to normalize.");

  if (is.null(bias)) {
  } else if (is.numeric(bias)) {
    bias <- rep(bias, length.out=length(fields));
  } else {
    throw("Invalid value of argument 'bias': ", bias);
  }
  bias0 <- bias;
  
  res <- list();
  for (slide in slides) {
    xy <- getSpotPosition(this, slide=slide);
    if (is.null(xy)) {
      layout <- getLayout(this);
      xy <- getPosition(layout);
      x <- xy[,1];
      y <- xy[,2];
    } else {
      x <- xy$x;
      y <- xy$y;
    }
    x0 <- x;
    y0 <- y;
    
    fieldIdx <- 0;
    for (field in fields) {
      fieldIdx <- fieldIdx + 1;
      z <- this[[field]][,slide];
      ok <- !is.na(z) & !is.infinite(z);
      x <- x0[ok];
      y <- y0[ok];
      z <- z[ok];
      if (length(z) == 0) {
        warning("No (valid) data points are available for field ", field, " on slide ", slide, ".");
        next;
      }
    
      cat("Spatial normalization of field ", field, " on slide ", slide, ". Please wait...", sep="");
      time <- system.time(
        if (method == "loess") {
          fit <- loess(z ~ cbind(x,y), family="symmetric", span=smooth, ...);
          zres <- as.vector(residuals(fit));
          zfit <- fitted(fit);
          rm(fit);
          bias <- bias0;
          if (is.null(bias[fieldIdx]))
            bias[fieldIdx] <- mean(zfit, na.rm=TRUE);
#          str(zres);
#          str(bias[fieldIdx]);
          this[[field]][ok, slide] <- zres + bias[fieldIdx];
          rm(zres);
          fitted <- rep(NA, nbrOfSpots(this));
          fitted[ok] <- zfit;
          res[[field]] <- cbind(res[[field]], fitted);
          dimnames(res[[field]]) <- NULL;
          rm(zfit);
        } else if (method == "ssanova") {
        }
      ); # time <- system.time(
      cat("done! (", formatC(time[3], format="f", digits=2), " seconds)\n", sep="");
    } # for (field in fields)
  } # for (slide in slides)

  clearCache(this); 

  invisible(res);
}, private=TRUE, trial=TRUE)




#########################################################################/**
# @set "class=MicroarrayData"
# @RdocMethod normalizePlatewise
#
# @title "Normalization performed plate by plate"
#
# \description{
#   Performs a normalization plate by plate on each of the slide seperately.
#   For details on plate-wise normalization see [1] and for details
#   on intensity dependent normalization see [2].\cr
#   
#   \emph{Note that
#   the data in the object is replaced with the new normalized data and
#   the old data is removed}. To keep the old data, make a copy of the
#   object before normalizing by using \code{clone(ma)}, see 
#   @see "R.oo::clone.Object" and example below.\cr
# }
#
# @synopsis
#
# \arguments{
#  \item{field}{The data field to be normalized.}
#  \item{method}{The normalization method to be used. If \code{"constant"}
#   the values of the data field will be shifted to have zero bias within
#   each plate group. If \code{"A"}, the values within each plate group
#   will be normalized against intensity dependent effects.}
#  \item{slides}{Slides to be included in the normalization. If @NULL,
#     all slides are normalized.}
#  \item{...}{Other arguments accepted by underlying normalization methods.}
# }
#
# \details{
#  The data will be normalize within each plate group \emph{individually}
#  by either a) assuming same shift and scale for whole plate i) shifting 
#  the data to have median/mean zero. ii) rescale the data so all the
#  groups have the same median absolute deviation (MAD) or standard
#  deviation (sd) or by b) assuming intensity dependent effects. For
#  more details see [1].
# }
#
# \examples{
#   # The option 'dataset' is used to annotate plots.
#   options(dataset="sma:MouseArray")
#
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   # Scaled intensity normalization print-tip by print-tip
#   ma.norm1 <- clone(ma)
#   normalizeWithinSlide(ma.norm1, method="s")
#
#   # Intensity normalization plate by plate
#   ma.norm2 <- clone(ma.norm1)
#   normalizePlatewise(ma.norm2, field="M", method="A")
#
#   # Plot data before and after normalization.
#   layout(matrix(1:9, ncol=3, byrow=TRUE))
#   plot(ma)
#   plotSpatial(ma)
#   plotPrintorder(ma)
#   plot(ma.norm1)
#   plotSpatial(ma.norm1)
#   plotPrintorder(ma.norm1)
#   plot(ma.norm2)
#   plotSpatial(ma.norm2)
#   plotPrintorder(ma.norm2)
# }
#
# \references{
#  \item{[1]}{Henrik Bengtsson, Plate Effects in cDNA microarray data, 
#    Matemathical Statistics, Centre for Matematical Sciences,
#    Lund University, Sweden. Manuscript, 2002.}
#  \item{[2]}{S. Dudoit, Y. H. Yang, M. J. Callow, and T. P. Speed. Statistical
#    methods for identifying differentially expressed genes in
#    replicated cDNA microarray experiments (Statistics, UC Berkeley,
#    Tech Report 578). 
#    URL: \url{http://www.stat.berkeley.edu/users/terry/zarray/Html/papersindex.html}}
# }
#
# @author
#
# \seealso{
#   For other within-slide normalization see @see "MAData.normalizeWithinSlide".
#   For across-slide normalization see @see "MAData.normalizeAcrossSlides".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("normalizePlatewise", "MicroarrayData", function(this, field, method="constant", slides=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);

  # Get information about which spots belongs to which plates
  plates <- getPlateGroups(getLayout(this))
  spots <- getSpots(plates);
  ngroups <- nbrOfGroups(plates);

  # Get the data to be normalized
  data <- this[[field]];

  # And if it is an intensity dependent normalization...
  if (method == "A") {
    if (!is.element("A", getFieldNames(this)))
      throw("Can not do intensity dependent normalization since this object does not have a field named 'A': ", data.class(this));
    if (field == "A")
      throw("It does not make sence to do intensity dependent normalization on field 'A' (intensity).");
    A <- this[["A"]];
  }

  # Do the normalization slide by slide
  for (slide in slides) {
#    cat("Normalizing slide ", slide, "...\n", sep="");
    y <- data[,slide];
    if (method == "constant") {
      y <- MicroarrayData$normalizeGroupsConstant(y, spots, ...)
    } else if (method == "A") {
      x <- A[,slide];
      y <- MicroarrayData$normalizeGroupsShiftFunction(x, y, spots, ...);
    } else {
      throw("Unknown plate-wise normalization method: ", method);
    }
    # Update the temporary data
    data[,slide] <- y;
  } # for (slide...)

  # Store the normalized data
  this[[field]] <- data;

  clearCache(this); 

  invisible(this);
})





setMethodS3("normalizeQuantile", "MicroarrayData", function(this, 
                                         fields, weights=1, robust=FALSE) {
  for (field in fields) {
    mat <- unclass(this[[field]]);
    mat <- normalizeQuantile(mat, weights=weights, robust=robust);
    idx <- 1:length(mat);
    this[[field]][idx] <- mat;
  }

  clearCache(this); 
}) # normalizeQuantile()




############################################################################
# HISTORY:
# 2005-03-23
# o Updated normalizeGroupsShiftFunction() so that approx() does not give 
#   warnings about 'Collapsing to unique x values' when doing lowess 
#   normalization.
# 2002-11-13
# o Updated normalizePlatewise() with the argument 'slides'.
# o Updated normalizeQuantile() to work faster for the new MicroarrayArray
#   class.
# 2002-10-24
# o Added normalizeQuantile().
# 2002-09-30
# o Added Rdoc comments for normalizePlatewise() and made it public.
# 2002-05-06
# * Added normalizePlatewise().
# * Implemented the two help functions for normalization of "any kind":
#     i) normalizeGroupsConstant()
#    ii) normalizeGroupsShiftFunction()
# * Extracted all normalization methods into MicroarrayData.NORM.R
# * Added a assertion that argument layout in the constructor is of class
#   Layout.
# 2002-05-04
# * Added applyGroupwise() and the applyXXXwise()'s.
# 2002-05-03
# * Added getBlanks().
# 2002-04-28
# * BUG FIX: adjustBiasScale...() did not scale correctly if bias was NOT
#   adjusted!
# 2002-04-21
# * Extracted log functions to MicroarrayData.LOG.R.
# * Extracted I/O functions to MicroarrayData.IO.R.
# * Extracted plot functions to MicroarrayData.PLOT.R.
# * Made getColors() generate grayscale colors by default. Before the method
#   was declared abstract.
# * Added trial version of normalizeGenewise().
# * Replaced some of the throw()'s with throw()'s.
# * Added getGeneReplicateIndex() and getGeneSlideReplicateIndex()
#   to MicroarrayData for fast access to the (gene, slide, replicate)
#   indices.
# 2002-04-20
# * Added trial versions of the static functions dataFrameToList() and
#   readToList(). These can support the read()  functions in the subclasses.
# * hist() now also excludes Inf's in addition to NA's. It turned out that
#   Lei Jiang's data, who reported the bug,  contained *one* M value that
#   was Inf. The result was that pretty(), which is called by hist(), gave
#   an error like "NA/NaN/Inf in foreign function call (arg 2)".
# * Added reference to 'plot' also in addition to 'par' in plot-functions.
#   This was done on a question how to set the limits on the axis.
# 2002-04-12
# * Updated the Rdoc example for plot() so it gives example on how to plot
#   a certain slide.
# 2002-04-06
# * Added support for multiple fields in normalizePrintorder() and
#   normalizeSpatial().
# * Renamed argument 'what' to 'field' in plotPrintorder().
# * Added argument 'breakpoints' to plotPrintorder().
# * Added normalizePrintorder().
# 2002-04-05
# * Added normalizeSpatial(). It is "smart", because it tries to get the
#   physical positions of the spots, but such information is not available
#   the positions according to the Layout object is used.
# 2002-04-03
# * write() now supports File objects.
# 2002-03-29
# * Updated the Rdoc's so any references to old get() are now to extract().
# 2002-03-25
# * Added static method createColors().
# 2002-03-24
# * Updated the Rdoc example for text().
# * BUG FIX: text() did not work correctly if argument 'labels' where not
#   explicitly set. The reason for this was that getInclude() was changed
#   to return a vector of TRUE and FALSE. Had to use which().
# * Changed this$getInclude(...) to getInclude(this, ...).
# 2002-03-13
# * BUG FIX: Forgot the sep="\t" in write().
# * Added some more Rdoc comments to write().
# 2002-03-10
# * BUG FIX: Argument 'slides' in write() was mistakenly named 'slide'.
# * write() in MicorarrayData will now by default save the object as the
#   data frame returned by as.data.frame(), i.e. if write() is not
#   overridden by a subclass e.g. GenePixData, which then saves in a
#   different format.
# * BUG FIX: extract() in MicroarrayData would neglect the special fields
#   "slide", "spot" and "gene". This automatically also fixed the fact that
#   as.data.frame() would not create these fields.
# 2002-02-26
# * Added the virtual fields "slide", "spot", "gene" to extract.
# 2002-02-25
# * Added read(), readAll(), write(), append().
# * Modified this class to support GenePixData etc directly without making
#   use of a ResultsData class.
# 2002-01-24
# * Renamed method get() to extract(). get() is not safe!
# * Rewritten to make use of setClassS3 and setMethodS3.
# 2002-01-19
# * Added getSlideName() and setSlideName(). Used in automatic labelling of
#   plots. See the putSlide() method.
# 2002-01-17
# * Added argument 'new=TRUE' to printReplicates. With new==FALSE it is
#   possible to plot another sequence of replicates in the same plot. This
#   is useful for instance when you want to look at the effect before and
#   after normalization.
# * BUG FIX: putGene() crashed if 'id' or 'name' was "auto"; instead of 
#   doing "if(id && name) ..." it is safer/better to do 
#   "if (id == TRUE && name == TRUE) ...".
# * BUG FIX: plotReplicates() didn't work if there where no within-slide
#   replicates. This is now corrected.
# 2002-01-15
# * Added putSlide().
# * Added seq() to simplify multiplots.
# * Added argument 'adjustMargins=TRUE' to subplots().
# * Added putTimestamp, putDataset, and putAuthor to all plot 
#   methods.
# 2001-11-18
# * Added getField() and getFields().
# 2001-08-11
# * Added plotPrintorder.
# 2001-08-08
# * Added the first core functionality of has-, get- setExcludedSpots.
#   By using which and unwhich I hope it could also be memory efficient.
#   I do not want to just set the values to NA to exclude, because then
#   one can not unexclude. Instead I want to flag some spots to be excluded.
# 2001-08-06
# * BUG FIX: plotSpatial set the plot history to "plotXY" instead of
#   "plotSpatial". This bug was probably a cut-and-paste misstake.
# * BUG FIX: plotXY and plotSpatial generated the wrong colors for all
#   cases where slide > 1. It turned out do be bug in getColors.MAData.
# 2001-08-03
# * Added getFieldNames(), setFieldNames(), renameField().
# 2001-07-31
# * Added hasLayout() and made setLayout() assert the argument layout.
# * Updated nbrOfSpots() to either ask the Layout of use 
#   get(fields=1, slides=1) to figure out the number of spots. This
#   method is a little bit inefficient so it should be overloaded by
#   subclasses.
# 2001-07-24
# * lowessCurve() now returns the lowess line.
# 2001-07-18
# * Bug fix: Forgot the 'labels' argument in call to text() in plotXY.
# 2001-07-17
# * TODO: Have to decide if spot indices should also run across slides, i.e.
#   the unique indices should continue counting on the next slide etc. This
#   is a complicated issue and somehow the though has to be digested before
#   making a decision. Now, I think some function are a little bit ambigous
#   in the use of slide, include and exclude arguments. It is pretty obvious
#   though that the include and exclude should be done before applying the
#   slide argument.
# * Added some Rdoc comments.
# * Updated the plot() method internally; now less if and then statements.
# * Removed the .lastPlot field. Making use of Device$setPlotParameters
#   instead.
# 2001-07-15
# * getIncludes() now also accepts lists in include/exclude arguments.
# * From now on the cex, col, pch etc arguments to the plot functions are
#   not subject to include and exclude. I made this decision since most
#   often you want for instance highlight four spots with four different
#   colors and nothing else. This was much harder to do before.
# 2001-07-12
# * Updated the pin.lty in plotXY to be the same as the one in Layout$put().
# * Bug fix: Trying to load a gpr data set with layout 8x4x15x16, the
#   pin.lty <- rep(...) function gave an error. used ngrid.c instead of
#   ngrid.r. The data I've tried this far have had ngrid.c == ngrid.r!
# 2001-07-11
# * Made .layout public, i.e. renamed it to layout.
# * Updated some of the Rdoc comments.
# 2001-07-09
# * Totally removed the use of image() in the plotSpatial() method. image()
#   had a "uncontrolable" color method and thanks to a highly improved
#   getPosition.Layout speed.
# 2001-07-06
# * Added addFlag(). Updated clearFlag() to work with regexpr's too.
# * Renamed plotYvsX to plotXY.
# 2001-07-05
# * Made include in the plot functions be operating on flags which have been
#   set by flag(). Removed all exlude from the plot functions.
# * Made plotYvsX() and plotSpatial() more generic and moved both of them to
#   this class. Also moved highlight(), text(), plot() to this class.
# 2001-07-04
# * Added the .extra field w/ the setExtra(), getExtra() methods.
# 2001-07-01
# * Removed plotSpatial(); now MicroarrayData is totally plot free.
# * Added getLabel() and setLabel(). Really useful!
# * Removed the rename() method since the internal field names should never
#   be modified.
# * Generic get() and set() seems to works fine. Added abstract setField().
# 2001-06-29
# * Created from old ResultsData.
############################################################################
