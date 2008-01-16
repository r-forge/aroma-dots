############################################################################
############################################################################
## 
##  NORMALIZATION METHODS
## 
############################################################################
############################################################################




#########################################################################/**
# @set "class=MAData"
# @RdocMethod normalizeAffine
#
# @title "Affine normalization based on non-logged data"
#
# \description{
#   @get "title".
#   For details, see @see "normalizeAffine.RGData".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "normalizeAffine.RGData".}
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
#   # Clone the data to get one non-normalized and one normalized data set.
#   maCurveFit <- clone(ma)
#   maAffine <- clone(ma)
#
#   # Normalize the data within slides using scaled print-tip normalization.
#   normalizeCurveFit(maCurveFit, groupBy="printtip", method="lowess")
#   normalizeAffine(maAffine, groupBy="printtip")
#
#   # Plot data before and after normalization.
#   subplots(9, nrow=3)
#   for (obj in list(ma, maCurveFit, maAffine)) {
#     # Plot M vs A and M spatially for array 1.
#     plot(obj)
#     plotSpatial(obj)
#     # Plot the densities of A for *all* arrays.
#     plotDensity(obj, what="A", xlim=c(4,16))
#   }
# }
#
# \seealso{
#   @seeclass
# }
#
# \author{
#   @get "author".
# }
#*/#########################################################################
setMethodS3("normalizeAffine", "MAData", function(this, ...) {
  rg <- as.RGData(this);
  normalizeAffine(rg, ...);
  ma <- as.MAData(rg);
  rm(rg);
  this$A <- ma$A;
  this$M <- ma$M;
  rm(ma);
  invisible(this);
}) # normalizeAffine()





#########################################################################/**
# @set "class=MAData"
# @RdocMethod normalizeWithinSlide
#
# @title "Within-slide normalization"
#
# \description{
#   Performs a within-slide normalization slide by slide. 
#   For a detailed explanation of normalization, see [1].\cr
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
#  \item{method}{The normalization method to be used. Currently there
#   are four different methods; 
#   \code{"m"} - median normalization, which sets the median of log
#   intensity ratios to zero, 
#   \code{"l"} - global lowess normalization,
#   \code{"p"} - print-tip group lowess normalization, and
#   \code{"s"} - scaled print-tip group lowess normalization.}
#  \item{weights}{Weights between zero and one, that is, in [0,1], of each
#   data point specifying how much that data points will affect the 
#   normalization. A data point with weight zero (or @FALSE) will not affect
#   the normalization, \emph{but} will be normalized. Currently only 0-1
#   (or @FALSE-@TRUE) weights are supported. Non-zero weights are treated
#   as ones.}
#  \item{lowess}{When doing global lowess normalization, \code{method="l"},
#   it is possible to specify the lowess line to be used. It is possible to
#   specify individual lowess lines for each of the slides by letting
#   lowess be a \emph{list} of lines. If \code{lowess=NULL} the lowess
#   curve will be estimated from the data.}
#
#   Note that \emph{only} one normalization is needed, i.e. doing different 
#   normalizations in serie on the same data set will not affect the
#   results. 
#
#   Also note that it is only the log ratios, M, are affected by the
#   normalization, i.e. the log intensities, A, are \emph{not} changed.
# }
#
# \references{
#   \item{[1]}{S. Dudoit, Y. H. Yang, M. J. Callow, and T. P. Speed. Statistical
#    methods for identifying differentially expressed genes in
#    replicated cDNA microarray experiments (Statistics, UC Berkeley,
#    Tech Report 578). 
#    URL: \url{http://www.stat.berkeley.edu/users/terry/zarray/Html/papersindex.html}}
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
#   # Clone the data to get one non-normalized and one normalized data set.
#   ma.norm <- clone(ma)
#
#   # Normalize the data within slides using scaled print-tip normalization.
#   normalizeWithinSlide(ma.norm, "s")
#
#   # Plot data before and after normalization.
#   layout(matrix(1:4, ncol=2, byrow=TRUE))
#   plot(ma)
#   plotSpatial(ma)
#   plot(ma.norm)
#   plotSpatial(ma.norm)
# }
#
# \note{
#   Note that the layout must be set for print-tip (\code{method="p"}) and
#   scaled (\code{method="s"}) normalization. If layout is not set, an
#   exception will be thrown. Normally, the layout is already set, such as
#   when the data is read from for instance GenePix, ScanAlyze and Spot.
# }
#
# \seealso{
#   For across-slide normalization see @seemethod "normalizeAcrossSlides".
#   @seeclass
# }
#
# \author{
#   @get "author".
#   Initial code for support of 'weights' by
#   Jon McAuliffe, Statistics Dept, UC Berkeley.
#   The original code was written by the sma authors
#    Yee Hwa Yang \email{yeehwa@stat.berkeley.edu}
#    Sandrine Dudoit \email{sandrine@stat.berkeley.edu} and
#    Natalie Roberts \email{nroberts@wehi.edu.au}.
# }
#*/#########################################################################
setMethodS3("normalizeWithinSlide", "MAData", function(this, method, slides=NULL, weights=NULL, lowess=NULL, f=0.3, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  slides <- validateArgumentSlides(this, slides=slides);

  # Argument 'method':
  if (missing(method))
    throw("Argument 'method' is missing.");

  methods <- rep(method, length.out=nbrOfSlides(this));
  for (method in methods) {
    isUnknown <- !is.element(method, c("l", "m", "p", "s", "printorder"));
    if (isUnknown)
      throw("Unknown normalization method: ", method);
  }

  # If the weights are TRUE/FALSE convert them into 1.0/0.0.
  if (is.null(weights)) {
    include <- TRUE;
    weights <- 1;
  } else if (is.numeric(weights)) {
    # Assert that the weights are valid.
    weights <- validateArgumentWeights(this, weights=weights);
    if (any(weights > 0.0 & weights < 1.0)) {
      warning("Currently only weights {0,1} are supported. Non-zero weights are set to one.");
    }
    include <- as.logical(weights);
  } else if (is.logical(weights)) {
    include <- weights;
    weights <- validateArgumentWeights(this, weights=weights);
  } else {
    throw("Unsupported class on argument 'weights': ", data.class(weights));
  }
  
  # Verify that argument 'lowess' is a structure with fields 'x' and 'y'
  if (length(lowess) > 0) {
    if ( any(!is.element(c("x","y"), names(lowess))) ) {
      throw("Can not perform global lowess normaliziation. Argument 'lowess' does not contain the field 'x' or 'y'.");
    } else if (any(!is.element(methods, c("l", "printorder")))) {
      throw("When argument 'lowess' is specified, the method must be 'l'.");
    }
  }

  # Asserts that all needed parameters are specified...
  layout <- getLayout(this);
  if (any(is.element(methods, c("p", "s", "printorder"))) && is.null(layout))
    throw("Can not normalize. Layout object is missing.");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For normalizing with a given lowess line
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  lowessNormalization <- function(A, M, include=TRUE, line=NULL, ...) {
    isFinite <- (is.finite(A) & is.finite(M));
    if (is.null(line)) {
      # Estimate lowess curve
      line <- lowess(A[isFinite & include], M[isFinite & include], ...);
    }
    # Adjusts only finite values - all other remains
    M[isFinite] <- M[isFinite] - approx(line, xout=A[isFinite], ties=mean)$y;
    M
  }
  
  printTipNormalization <- function(A, M, include=TRUE, layout, ...) {
    npin <- nbrOfGrids(layout);
    pin <- c(0, rep(gridSize(layout), npin) * (1:npin));
    include <- rep(include, length.out=length(M));
    idx <- 1:length(M);
    for(j in 1:npin) {
      index <- ((pin[j] + 1) <= idx) & (idx <= pin[j + 1])
      M[index] <- lowessNormalization(A[index], M[index], include=include[index], ...)
    }
    M
  }

  printTipScaleNormalization <- function(M, layout) {
    # 1. Arange the data in a matrix where data from the same printtip is
    #    in the same column. M: [1:gridSize, ]
    M <- matrix(M, nrow=gridSize(layout));
    
    # 2. Get the mad for each print-tip group. M.mad: [1:gridSize]
    M.mad <- apply(M, MARGIN=2, FUN=mad, na.rm=TRUE);

    # 3. Scale each print-tip group with a print-tip weight.
    # 3a. Average M.mad
    l.M.mad <- log(M.mad);
    sum.M.mad <- sum(l.M.mad[!is.infinite(l.M.mad)], na.rm=TRUE)
    
    # 3b. Print-tip weights
    weight <- M.mad / exp( sum.M.mad / nbrOfGrids(layout) );

    # 3c. Scale each print-tip group with the print-tip weights
    scaledM <- t(t(M) / weight);

    # 4. Rearrange data to original shape
    matrix(scaledM, nrow=nbrOfSpots(layout))
  }

 
  K <- length(lowess);
  if (K == 1 && all(is.element(c("x","y"), names(lowess))))
    lowess <- list(lowess);

  # Normalize slide by slide.
  for (kk in slides) {
    # Get the current (A,M) data
    M <- this$M[,kk];
    A <- this$A[,kk];

    if (methods[kk] == "m") {
      # Normalize so the median of the M-vales becomes zero.
      M <- M - median(M, na.rm=TRUE);
      description <- "median"
    } else if (methods[kk] == "l") {
      line <- if (length(lowess) > 0) lowess[[(kk-1) %% K + 1]] else NULL;
      M <- lowessNormalization(A, M, include=include, line=line, f=f, ...);
      description <- "lowessed";
    } else if (methods[kk] == "p" || methods[kk] == "printtip") {
      M <- printTipNormalization(A, M, layout, include=include, f=f, ...);
      description <- "print-tip"
    } else if (methods[kk] == "s") {
      M <- printTipNormalization(A, M, layout, include=include, f=f, ...);
      M <- printTipScaleNormalization(M, layout);
      description <- "scaled print-tip"
    } else if (methods[kk] == "printorder") {
      tidx <- toPrintorderMatrix(layout);
      description <- "printorder";
      # Order the data in print order
      Mt <- M[tidx];
#      At <- A[tidx];
      idx <- 1:nbrOfSpots(this);
      # Now, do regular lowess normalization
#      ind <- !is.finite(At);
      isFinite <- is.finite(Mt);
      if (is.null(lowess)) {
        includet <- include[tidx];
        lowess <- lowess(idx[isFinite & includet], Mt[isFinite & includet], ...);
        rm(includet);
      }
      Mt[isFinite] <- Mt[isFinite] - approx(lowess, xout=idx[isFinite], ties=mean)$y;
#      At[isFinite] <- At[isFinite] - approx(lowess, xout=idx[isFinite], ties=mean)$y;
      # Update the normalized values
      M[tidx] <- Mt;
#      A[tidx] <- At;
#      this$A[,kk] <- A;
    }
    this$M[,kk] <- M;
  } # for (kk in 1:nbrOfSlides(this))

  # If labels are set, update the with a description.
  label <- getLabel(this, "M");
  if (!is.null(label)) {
    if (is.expression(label))
      label <- label[[1]];
    label <- substitute(paste(X, ", ", plain(Y)), list(X=label, Y=description));
    setLabel(this, "M", label);
  }

  # Data object has changes, clear cached values.
  clearCache(this);
  
  invisible(this);
})


#########################################################################/**
# @RdocMethod normalizeAcrossSlides
#
# @title "Normalizes across slides"
#
# @synopsis
#
# \arguments{
#   \item{slides}{The set slides to be used and to be normalized. If
#     @NULL all slides will be normalized.}
#   \item{newMAD}{After the normalization all slides will have an maximum
#     absolute deviation of \code{newMAD}. If @NULL, the result will
#     be the same as if \code{newMAD} would be set to the geometrical mean
#     of each individual slide's MAD.}
# }
#
# \description{
#   Normalizes across some or all slides. After doing within-slide normalization
#   (see @seemethod "normalizeWithinSlide", an across-slide normalization must 
#   be performed before the data on the different slides can be compared with
#   each. Across-slide normalization scales the log ratios (M) for all slide so
#   each slide gets the same log-ratio deviation based on the robust deviation
#   measure Maximum Absolute Deviation (MAD).
#
#   If one would like to normalize the deviation across \emph{different data set},
#   i.e. different @see "MAData" objects, one can make use of the argument 
#   \code{newMAD}, which forces the slides to get a specific median absolute
#   deviation value.
#
#   Note, in the case where one set of slides comes from one type of experimental
#   setup and a second set of slides comes from another setup, and they are stored
#   in the same \code{MAData} object, these two groups of slides \emph{can} be 
#   normalized together using one (in other words, you do \emph{not} have to
#   normalize the two groups seperately and the rescale them with \code{newMAD}).
#
#   Also note that it is only the log ratios, M, are affected by the
#   normalization, i.e. the log intensities, A, are \emph{not} changed.
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
#   layout(matrix(1:9, ncol=3, byrow=TRUE))
#
#   # Plot data before within-slide normalization
#   for (k in 1:3)
#     boxplot(ma, groupBy="printtip", slide=k, main=paste("No normalization - #", k, sep=""))
#
#   # Plot data after scaled print-tip normalization
#   normalizeWithinSlide(ma, "s")
#   for (k in 1:3)
#     boxplot(ma, groupBy="printtip", slide=k, main=paste("Within-slide norm. - #", k, sep=""))
#
#   # Plot data after across-slide normalization
#   normalizeAcrossSlides(ma)
#   for (k in 1:3)
#     boxplot(ma, groupBy="printtip", slide=k, main=paste("Across-slide norm. - #", k, sep=""))
# }
#
# @author
#
# \seealso{
#   For an detailed explanation of the robust deviation measure MAD 
#   (Median Absolute Deviation) see @see "stats::mad".
#   For within-slide normalization see @seemethod "normalizeWithinSlide".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("normalizeAcrossSlides", "MAData", function(this, slides=NULL, newMAD=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  slides <- validateArgumentSlides(this, slides=slides);

#  if (length(slides) == 1) return(invisible()); # Nothing will be changed!
  
  # Calculate Median Absolute Deviation (MAD) for each slide.
  M <- as.matrix(this$M[, slides]);
  M.MAD <- apply(M, MARGIN=2, FUN=mad, na.rm=TRUE);
  M2.MAD <- M.MAD[is.finite(M.MAD)];
  n <- length(M2.MAD);
  if (n != length(slides))
    warning("Some of the slides have a MAD value which is NA, NaN or Inf, which is really strange!");

#  if (n == 1) return(invisible()); # Nothing will be changed!

  # Either use 'newMAD' or use the geometrical mean of MAD's...
  if (!is.null(newMAD)) {
    Z <- newMAD;
  } else {
  	# Calculate the all-slide normalization factor Z as the geometrical mean
  	# of all slides' MADs.
  	Z <- prod(M2.MAD)^(1/n);  # (a scalar)
  }

  # Calculate the scale factor for each slide.
  scale <- Z / M.MAD;

  # Scale each slide individually, by its scale factor.
  for (k in seq(slides)) {
    this$M[,slides[k]] <- this$M[,slides[k]] * scale[k];
  }

  label <- getLabel(this, "M");
  if (is.null(label))
    label <- expression(M==log[2](R/G));
  if (is.expression(label))
    label <- label[[1]];
  description <- "across";
  label <- substitute(paste(X, ", ", plain(Y)), list(X=label, Y=description));
  setLabel(this, "M", label);

  # Data object has changes, clear cached values.
  clearCache(this);

  invisible(this);
})



setMethodS3("normalizeGenewise", "MAData", function(this, fields="M", bias=0, scale=1, ...) {
  normalizeGenewise.MicroarrayData(this, fields=fields, bias=bias, scale=scale, ...);
})



setMethodS3("normalizePrintorder", "MAData", function(this, fields=c("M","A"), ...) {
  normalizePrintorder.MicroarrayData(this, fields=fields, ...);
})



setMethodS3("normalizePlatewise", "MAData", function(this, field="M", ...) {
  normalizePlatewise.MicroarrayData(this, field=field, ...);
})








setMethodS3("normalizeSpatial", "MAData", function(this, fields=c("M","A"), ...) {
  normalizeSpatial.MicroarrayData(this, fields=fields, ...);
})




setMethodS3("normalizeQuantile", "MAData", function(this, 
                                                  fields=c("M", "A"), ...) {
  normalizeQuantile(this, fields=fields, ...);
}) # normalizeQuantile()




#########################################################################/**
# @RdocMethod normalizeLogRatioShift
#
# @title "Within-slide normalization that adjust the log-ratios shift"
#
# @synopsis
#
# \arguments{
#   \item{slides}{The slides which should be included in the calculations.
#     If @NULL, all slides are included.}
#   \item{groupBy}{@character string or @see "LayoutGroups" specifying the
#     groups of spots that are to normalized individually. 
#     If @NULL, all spots are normalized together.}
#   \item{method}{@character string specifying if the median or the mean
#     should be used to estimate the shift of the log-ratios.}
# }
#
# \description{
#   @get "title".
# }
#
# \value{
#   Returns a @list structure containing information about the fit etc.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("normalizeLogRatioShift", "MAData", function(this, slides=NULL, groupBy=NULL, method=c("median", "mean")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  slides <- validateArgumentSlides(this, slides=slides);
  spots <- validateArgumentGroupBy(this, groupBy=groupBy);
  method = match.arg(method);
  if (method == "median") {
    fcn <- function(x) median(x, na.rm=TRUE);
  } else if (method == "mean") {
    fcn <- function(x) mean(x, na.rm=TRUE);
  }

  fit <- NULL;
  as <- matrix(NA, nrow=nbrOfSpots(this), ncol=nbrOfSlides(this));
  for (slide in slides) {
    Mall <- this$M[,slide];

    for (spot in spots) {
      M <- Mall[spot];
      a <- fcn(M);
      Mall[spot] <- M - a;

      # Store the amount of shift for each spot
      as[spot,slide] <- a;
    } # for (spot ...)

    this$M[,slide] <- Mall;
  } # for (slide ...)

  fit <- list(shift=as);
  invisible(fit);
}, protected=TRUE); 




#########################################################################/**
# @RdocMethod estimateMACurve
#
# @title "Estimates a smooth intensity-dependent curve in (A,M)"
#
# @synopsis
#
# \arguments{
#   \item{slides}{The slides which should be included in the calculations.
#     If @NULL, all slides are included.}
#   \item{groupBy}{@character string or @see "LayoutGroups" specifying the
#     groups of spots that are to normalized individually. 
#     If @NULL, all spots are normalized together.}
#   \item{weights}{A @vector or @matrix of spot weights used when 
#     estimating the normalization function.}
#   \item{method}{@character string specifying which method to use when
#     fitting the intensity-dependent function.
#     Supported methods: 
#      \code{"loess"} (better than lowess), 
#      \code{"lowess"} (classic; supports only zero-one weights),
#      \code{"spline"} (more robust than lowess at lower and upper
#                       intensities; supports only zero-one weights),
#      \code{"robust.spline"} (better than spline).
#   }
#   \item{span}{A @double value specifying the bandwidth of the estimator used.}
#   \item{...}{No used (allows other methods to pass additional garbage).}
# }
#
# \description{
#   @get "title".
# }
#
# \value{
#   Returns a @list where each element is a @list corresponding to a slide, 
#   which in turn contains estimates of each set of spots specified by
#   the \code{groupBy} argument.
# }
#
# \section{Missing values}{
#  The estimate will only be made based on complete finite observations.
# }
#
# @author
#
# \seealso{
#   Utilized by @seemethod "drawCurveFit" and @seemethod "normalizeCurveFit".
#
#   @seeclass
# }
#*/#########################################################################
setMethodS3("estimateMACurve", "MAData", function(this, slides=NULL, groupBy=NULL, weights=NULL, method=c("loess", "lowess", "spline", "robust.spline"), span=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  slides <- validateArgumentSlides(this, slides=slides);
  spots <- validateArgumentGroupBy(this, groupBy=groupBy);
  weights <- validateArgumentWeights(this, weights=weights);
  method = match.arg(method);

  if (method == "lowess") {
    weights <- validateArgumentWeights(this, weights=weights, zeroOneOnly=TRUE);

    fcn <- function(A, M, weights=NULL, span=0.3, ...) {
      incl <- if (!is.null(weights)) (weights > 0) else TRUE;
      fit <- doCall("lowess", x=A[incl], y=M[incl], f=span, ...);
      fit$predictM <- function(newA) approx(fit, xout=newA, ties=mean)$y;
      fit;
    }
    if (is.null(span))
      span <- 0.3;
  } else if (method == "loess") {
    fcn <- function(A, M, weights=NULL, span=0.75, ...) {
      fit <- doCall("loess", formula=M ~ A, weights=weights, 
                    family="symmetric", degree=1, span=span, 
                    control=loess.control(trace.hat="approximate", 
                                     iterations=5, surface="direct"), ...);
      fit$predictM <- function(newA) predict(fit, newdata=newA);
      fit;
    }
    if (is.null(span))
      span <- 0.75;
  } else if (method == "spline") {
    weights <- validateArgumentWeights(this, weights=weights, zeroOneOnly=TRUE);

    fcn <- function(A, M, weights=NULL, span=0.75, ...) {
      incl <- if (!is.null(weights)) (weights > 0) else TRUE;
      fit <- doCall("smooth.spline", x=A[incl], y=M[incl], spar=span, ...);
      fit$predictM <- function(newA) predict(fit, x=newA)$y;
      fit;
    }
    if (is.null(span))
      span <- 0.75;
  } else if (method == "robust.spline") {
    fcn <- function(A, M, weights=NULL, span=0.75, ...) {
      require(R.basic) || throw("Could not load package R.basic.");
      fit <- doCall("robust.smooth.spline", x=A, y=M, w=weights, spar=span, ...);
      fit$predictM <- function(newA) predict(fit, x=newA)$y;
      fit;
    }
    if (is.null(span))
      span <- 0.75;
  }

  if (!is.numeric(span) || span <= 0 || span > 1)
    throw("Argument 'span' must be with [0,1): ", span);

  fit <- list();
  for (slide in slides) {
    Aall <- this$A[,slide];
    Mall <- this$M[,slide];
    if (is.matrix(weights)) {
      weightsAll <- weights[,slide];
    } else {
      weightsAll <- weights;
    }

    fitSlide <- list();
    for (kk in seq(along=spots)) {
      spot <- spots[[kk]];
      M <- Mall[spot];
      A <- Aall[spot];
      isFinite <- (is.finite(A) & is.finite(M));
      w <- if (is.null(weightsAll)) NULL else weightsAll[isFinite];
      fitSlide[[kk]] <- fcn(A[isFinite], M[isFinite], weights=w, span=span);
    } # for (spot ...)
    rm(A,M,w);

    fit[[slide]] <- fitSlide;
  } # for (slide ...)
  rm(Aall,Mall,weightsAll,fitSlide);

  names <- names(slides);

  invisible(fit);
}, protected=TRUE) # estimateMACurve()




#########################################################################/**
# @RdocMethod normalizeCurveFit
#
# @title "Within-slide normalization that adjust log-ratios by estimating a smooth intensity-dependent curve in (A,M)"
#
# @synopsis
#
# \arguments{
#   \item{slides}{The slides which should be included in the calculations.
#     If @NULL, all slides are included.}
#   \item{groupBy}{@character string or @see "LayoutGroups" specifying the
#     groups of spots that are to normalized individually. 
#     If @NULL, all spots are normalized together.}
#   \item{...}{Other arguments, such as \code{groupBy} and \code{method}, 
#      passed to the curve estimator @seemethod "estimateMACurve".}
# }
#
# \description{
#   @get "title".
# }
#
# \value{
#   Returns itself invisibly.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("normalizeCurveFit", "MAData", function(this, slides=NULL, groupBy=NULL, ..., .fits=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  slides <- validateArgumentSlides(this, slides=slides);
  spots <- validateArgumentGroupBy(this, groupBy=groupBy);

  fits <- .fits;
#  if (is.null(fits))
#    fits <- estimateMACurve(this, slides=slides, groupBy=groupBy, ...);
  rm(.fits);

  rg <- as.RGData(this);

  for (slide in slides) {
    fitSlide <- fits[[slide]];
    Xall <- cbind(rg$R[,slide], rg$G[,slide]);

    for (kk in seq(along=spots)) {
      spot <- spots[[kk]];
      fitSpot <- fitSlide[[kk]];
      X <- Xall[spot,];
      
      X <- normalizeCurveFit(X, ...);

      Xall[spot,] <- X;
    } # for (spot ...)
    rm(X);

    rg$R[,slide] <- Xall[,1];
    rg$G[,slide] <- Xall[,2];
  } # for (slide ...)
  rm(Xall);

  this$M <- rg$M;
  this$A <- rg$A;
print(rg);
  rm(rg);

  # Data object has changes, clear cached values.
  clearCache(this);

  invisible(this);
}) # normalizeCurveFit();







setMethodS3("drawCurveFit", "MAData", function(this, slides=NULL, groupBy=NULL, ..., .fits=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  slides <- validateArgumentSlides(this, slides=slides);
  spots <- validateArgumentGroupBy(this, groupBy=groupBy);

  fits <- .fits;
  if (is.null(fits))
    fits <- estimateMACurve(this, slides=slides, groupBy=groupBy, ...);
  rm(.fits);

  for (slide in slides) {
    fitSlide <- fits[[slide]];
    Aall <- this$A[,slide];

    for (kk in seq(along=spots)) {
      spot <- spots[[kk]];
      fitSpot <- fitSlide[[kk]];
      A <- Aall[spot];

      A <- A[is.finite(A)];
      newA <- seq(from=min(A), to=max(A), length=100);

      xy <- list(x=newA, y=fitSpot$predict(newA));
      lines(xy, ...);
    } # for (spot ...)
  } # for (slide ...)

  invisible(this);
}) # drawCurveFit()







#########################################################################/**
# @RdocMethod normalizeAffineShift
#
# @title "Within-slide normalization that adjust the channel bias"
#
# @synopsis
#
# \arguments{
#   \item{slides}{The slides which should be included in the calculations.
#     If @NULL, all slides are included.}
#   \item{groupBy}{@character string or @see "LayoutGroups" specifying the
#     groups of spots that are to normalized individually. 
#     If @NULL, all spots are normalized together.}
#   \item{method}{@character string specifying which method to use. 
#     If \code{"Kerr"}, perpendicular translation normalization is applied.
#     If \code{"Newton"}, parallel translation normalization is applied.
#     If \code{"Bengtsson"}, semi-perpendicular translation normalization
#     is applied.}
#   \item{interval}{The interval of possible shifts. 
#     Either a @numeric @vector of length 2 or a @character string.
#     If \code{"safe"}, the interval will be set, depending on method,
#     such that no signals are made non-positive.
#     If \code{"auto"}, the interval is choosen automatically depending on
#     method. For some methods, \code{"auto"} gives the same result
#     as \code{"safe"}.}
#   \item{...}{Currently not used.}
# }
#
# \description{
#   @get "title" with the objective to minimize intensity dependency of
#   the log-ratios. 
# }
#
# \value{
#   Returns a @list structure containing information about the fit etc.
# }
#
# @examples "../incl/MAData.normalizeAffineShift.Rex"
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("normalizeAffineShift", "MAData", function(this, slides=NULL, groupBy=NULL, weights=NULL, method=c("Kerr", "Newton", "Bengtsson"), interval=c("safe", "auto"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  slides <- validateArgumentSlides(this, slides=slides);
  spots <- validateArgumentGroupBy(this, groupBy=groupBy);
  weights <- validateArgumentWeights(this, weights=weights);
  method = match.arg(method);
  if (is.character(interval))
    interval = match.arg(interval);
  rg <- as.RGData(this);

  fit <- NULL;

  # sd  - non-robust, i.e. do not use this!
  # IQR - robust 
  # mad - more robust
  if (method == "Kerr") {
    # objective function to be minimized
    objectiveFunction = function(a, R,G) {
      M <- log((R+a)/(G-a), base=2);
      mad(M, na.rm=TRUE);
    }

    negBefore <- sum(rg$R < 0 | rg$G < 0, na.rm=TRUE);

    as <- matrix(NA, nrow=nbrOfSpots(this), ncol=nbrOfSlides(this));
    for (slide in slides) {
      Rall <- rg$R[,slide];
      Gall <- rg$G[,slide];

      for (spot in spots) {
        R <- Rall[spot];
        G <- Gall[spot];

  	# Get the search interval
  	if (is.character(interval)) {
  	  if (interval == "safe") {
  	    int <- c(-min(R,na.rm=TRUE),min(G,na.rm=TRUE))
  	  } else if (interval == "auto") {
  	    int <- c(-1,1)*max(R,G, na.rm=TRUE);
  	  }
  	} else {
  	  int <- interval;
  	}
  
  	# Find optimal shift 
  	a = optimize(objectiveFunction, interval=int, maximum=FALSE, R=R,G=G)$minimum;
  
  	# Count how many invalid data points the shift introduces.
  	R <- R + a;
  	G <- G - a;
  
  	# Apply optimal shift 
  	Rall[spot] <- R;
  	Gall[spot] <- G;

        # Store the amount of shift for each spot
        as[spot,slide] <- a;
      } # for (spot ...)

      rg$R[,slide] <- Rall;
      rg$G[,slide] <- Gall;
    } # for (slide ...)
    negAfter <- sum(rg$R < 0 | rg$G < 0, na.rm=TRUE);

    ma <- as.MAData(rg);
    rm(rg);
    this$M <- ma$M;
    this$A <- ma$A;
    rm(ma);

    newNegs <- (negAfter - negBefore);
    if (newNegs > 0)
      warning("Kerr's optimal shift normalization introduced negative signals.");
    fit <- list(shift=as, nbrOfLostValues=newNegs, objectiveFunction=objectiveFunction);
  } else if (method == "Bengtsson") {
    # objective function to be minimized
    objectiveFunction = function(a, R,G) {
      M <- log((R+a*(a>0))/(G-a*(a<0)), base=2);
      mad(M, na.rm=TRUE);
    }

    negBefore <- sum(rg$R < 0 | rg$G < 0, na.rm=TRUE);

    as <- matrix(NA, nrow=nbrOfSpots(this), ncol=nbrOfSlides(this));
    for (slide in slides) {
      Rall <- rg$R[,slide];
      Gall <- rg$G[,slide];

      for (spot in spots) {
        R <- Rall[spot];
        G <- Gall[spot];

  	# Get the search interval
  	if (is.character(interval)) {
  	  if (interval == "safe" || interval == "auto")
  	    int <- c(-1,1)*max(G,R,na.rm=TRUE);
  	} else {
  	  int <- interval;
  	}
  
  	# Find optimal shift 
  	a = optimize(objectiveFunction, interval=int, maximum=FALSE, R=R,G=G, tol=1e-2)$minimum;
  
  	# Count how many invalid data points the shift introduces.
  	# (This should never happen, but we check it anyway to assert
  	#  the correctness of the method).
  	R <- R + a * (a > 0);
  	G <- G - a * (a < 0);

  	# Apply optimal shift 
  	Rall[spot] <- R;
  	Gall[spot] <- G;

        # Store the amount of shift for each spot
        as[spot,slide] <- a;
      } # for (spot ...)

      rg$R[,slide] <- Rall;
      rg$G[,slide] <- Gall;
    } #for (slide ...)
    negAfter <- sum(rg$R < 0 | rg$G < 0, na.rm=TRUE);

    ma <- as.MAData(rg);
    rm(rg);
    this$M <- ma$M;
    this$A <- ma$A;
    rm(ma);

    newNegs <- (negAfter - negBefore);
    if (newNegs > 0) {
      throw(InternalErrorException("Bengtsson's optimal shift normalization introduced negative signals. Please, contact the author as this should never happen!", package=aroma));
    }
    fit <- list(shift=as, nbrOfLostValues=newNegs, objectiveFunction=objectiveFunction);
  } else if (method == "Newton") {
    throw("Newton's method is not implemented yet.");
  }

  # Data object has changes, clear cached values.
  clearCache(this);

  this$fit <- fit;

  invisible(fit);
}, private=TRUE);  # normalizeAffineShift()


############################################################################
# HISTORY:
# 2008-01-15
# o CLEAN UP: Removed non-used argument 'newA' from internal fit function
#   of estimateMACurve() in MAData.
# 2006-02-08
# o The code for normalizeWithinSlide() assumed 'method' was a single value,
#   but it may be a vector which then generates a warning.
# 2005-03-23
# o Updated normalizeWithinSlide() so that approx() does not give warnings
#   about 'Collapsing to unique x values' when doing lowess normalization.
# 2005-01-23
# o Now drawCurveFit() and normalizeCurveFit() utilizes estimateMACurve().
# 2004-12-28
# o Added normalizeAffine() with makes use ditto in RGData. Added an 
#   example which compares raw data, lowess, and affinely normalized data
#   in M vs A, spatially, and densities of A.
# o Renamed normalizeAffine() for MAData to normalizeAffineShift() to be
#   consistent with normalizeAffine() in RGData.
# o normalizeCurveFit() and drawCurveFit() updated and tested.
# o Updated the loess parameters for normalizeCurveFit() according to 
#   G. Smyth's suggestions by email on 2003-03-17. 
# 2004-12-27
# o Added validation of argument 'weights' for normalizeWithinSlide().
# 2004-01-12
# o Update code to make use of validateArgumentGroupBy().
# o Made normalizeAffine() for the MAData class private again. The one
#   in aroma.affine is much better.
# 2003-09-27
# o Added normalizeLogRatioIntensityBias() for lowess and loess curve-fit
#   normalization.
# o Renamed normalizeChannelBias() to normalizeAffine() to make it
#   more clear that it is a scalar bias, not an intensity dependent bias.
# o Created normalizeLogRatioShift() for median (or mean) log-ratio shifts.
# o Added groupBy argument to normalizeChannelBias(). Now we can normalize
#   by printtips, plates etc.
# 2003-09-18
# o Added protected normalizeChannelBias().
# 2003-04-08
# o Removed all "Default value is..." in the Rdoc comments.
# 2002-10-24
# o Added normalizeQuantile().
# 2002-09-30
# o Added wrapper method normalizePlatewise() with default field="M".
# 2002-06-30
# * Splitted MAData.R into MAData.R and MAData.NORM.R.
# * BUG FIX: getMOR() gave an error if probs==-0.5 ("mean") only.
# 2002-06-24
# * BUG FIX: mean() was broken.
# * Made getGeneDiscrepancies() obsolete.
# * Added argument 'weights' to normalizeWithinSlide() on request from Jon
#   McAuliffe. Currently only 0-1 weights or FALSE-TRUE weights are
#   supported, since lowess() does not support weights. If loess() would
#   be used weights could be in [0,1].
# 2002-05-26
# * Replaced all data( MouseArray) with SMA$loadData("mouse.data") because
#   it save memory.
# 2002-05-11
# * BUG FIX: When removed as.SSMatrix() forgot to replace with as.matrix().
#   Made as.TMAData() to fail for instance.
# 2002-05-10
# * Changed the default red-to-green range for getColors() from [-3,3] to
#   [-2,2]. This was done to create for colorful pictures. Before they were
#   a little bit too brownisch.
# 2002-05-07
# * Added plotDiporder()
# 2002-05-06
# * Updated getGeneDiscrepancies() with standardize=FALSE.
# 2002-05-05
# * Removed the transformation to class SSMatrix. Is it really needed? It
#   also made the data loose its slide name (colnames).
# * Added getGeneDiscrepancies(). This is a better measure of goodness than
#   getGeneResiduals() since in the end of the day we are only interested
#   in the genes and not the spots per se.
# * Added getGeneResiduals(). Idea is got get some kind of goodness measure
#   of the slide.
# 2002-04-21
# * Replaced a few throw()'s with throw()'s.           
# * Added normalizeGenewise() with bias=10 for the A's.
# 2002-04-20
# * Added a trial version of read(). Most of the job is done in support
#   functions in the MicroarrayData class.
# 2002-04-06
# * Added normalizeSpatial().
# * Added plot3d().
# * Added normalizePrintorder().
# * Added plotPrintorder().
# 2002-04-03
# * Updated normalizeAcrossSlides() to accept normalization of one single
#   slide. Some user has slides with different layout and this is the way
#   to normalize across such slides. They have to specify 'newMAD'.
# 2002-02-26
# * Removed extract(), as.data.frame(), nbrOfSlides() etc.
# * Updated the Rdoc's.
# 2002-01-24
# * Rewritten to make use of setClassS3 and setMethodS3.
# * Renamed all get() to extract().
# 2002-01-19
# * Added support for slide names in as.MAData(), as.RGData().
# * Added as.RawData() for being consistent with the RGData methods.
# 2002-01-18
# * Added argument 'newMAD' to force the scale of the slides to a certain
#   deviation.
# * normalizeWithinSlide now gives a descriptive error message if argument
#   'method' is missing.
# 2002-01-10
# * Internally, now M and A are SSMatrix's, which are easy to convert to
#   GSRArray's (indexed by Genes, Slides, Replicates).
# 2001-11-12
# * Added "MAData:" in as.character().
# * Typofix: The Rd example for as.RGData was never generated.
# 2001-09-30
# * Bug fix: If normalizeAcrossSlides() where called on an object with only
#   one slide, it gave an error. Now it return nicely instead with nothing
#   changed.
# 2001-09-29
# * Added as.BMAData (the Bayesian model of Lönnstedt et al.)
# * Added the argument 'slides' to normalizeAcrossSlides().
# 2001-09-28
# * Bug fix: Forgot the argument 'f' in normalizeWithinSlide(). This
#   resulted in a normalization that was not consistent with the line(s)
#   drawn by lowessCurve().
# 2001-08-10
# * Updated the normalizeWithinSlide so it does not depend on plot.mva in
#   sma anymore. Was this the last sma function to be replaced?
# 2001-08-08
# * Updates dyeSwap with the argument 'slides' and wrote its Rdoc comments.
# 2001-08-07
# * First version of having a optional 'lowess' argument to
#   normalizeWithinSlide.MAData.
# 2001-08-06
# BUG FIX: getColors() did not work with slide > 1.
# 2001-08-02
# * Updated mean() and var().
# 2001-08-01
# * Added the field spot in as.data.frame().
# 2001-07-14
# * Removed the default method in normalizeWithinSlide. As Javier mentioned,
#   having "n" (no normalization) by default as before, confuses people. 
#   This was an artifact from the plot.mva function in sma. By not, having
#   a default method, the user is forced to think about what normalization
#   method to use.
# * Update equals to include Layout comparison too. Added Rdoc for equals().
# 2001-07-11
# * Updated some of the Rdoc comments.
# * Bug fix: Mistakenly named normalizeAcrossSlides as normalized...().
# 2001-07-05
# * Removed plotMvsA(), plotMvsAGridwise(), and plotAbsMvsA().
# 2001-07-04
# * Removed ttest(). Now only as.TMAData() exists.
# 2001-07-02
# * Implemented the new support for highlighting.
# 2001-07-01
# * Removed the getNormalized*() methods. Better if user instead use clone.
# * Removed the .lastLayout field. Not needed anymore.
# * Updated the Rdoc comments alot.
# * Made a lot of the methods return "itself" if "itself" wasn't an 
#   reference. A little bit backward compatible with the standard [R]
#   object-oriented approach, which doesn't support reference.
# 2001-06-29
# * Bug fix: highlighting spots after plotSpatial flipped the x and the y
#   axis around. Now fixed.
# * Updated the Rdoc comments.
# * Added as.data.frame().
# 2001-06-08
# * Added las=3, i.e. style of axis labels is "always vertical", to 
#   plotBoxPlot() since "12", "14", "16" etc would normally not print.
# * Added getRange(). Useful for plotting several slide with the same scale.
# 2001-06-07
# * Added demo()
# * Bug fix: plotMvsA() and plotAbsMvsA() didn't pay attention to the
#   argument 'slide'.
# * Added arguments slide, include and exclude to getColors(). Improved
#   getColors() so it accepts NA's.
#   Argument col in plotMvsA() etc now supports palettes too, e.g. the value
#   "redgreen" will create red to green colored spots depending on M and A.
# 2001-06-04
# * Added as.matrix() to the constructor. 
# * Added getM(), setM(), getA() and setA().
# 2001-05-14
# * Added getInternalReferences() for improving gco() performance.
# 2001-04-06
# * Added normalizeAcrossSlides() and getNormalizedAcrossSlides().
# 2001-04-04
# * Added getDenseSpots().
# 2001-04-02
# * Implemented highlight() for plotSpatials too! However, it needs some
#   speedup improvements. Probably in layout class.
# * Added highlight() and made it "clever", i.e. it knows whats on the
#   x and the y axes.
# 2001-04-01
# * Added denominator argument to ttest.
# * Added support for stylers in all plot methods
# 2001-03-28
# * Updated normalize to normalize multiple slide at once.
# 2001-03-27
# * Added histogram() and getHistogram().
# * Question: Should normalize() creata a new MAData object or not?!?
#   Maybe add a method getNormalized() which does this?
# * Added normalize(), plotSpatial(), plotBox().
# 2001-03-10
# * Created.
############################################################################
