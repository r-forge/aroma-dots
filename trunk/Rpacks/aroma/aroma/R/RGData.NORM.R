#########################################################################/**
# @set "class=RGData"
# @RdocMethod normalizeAffine
#
# \encoding{latin1}
# 
# @title "Weighted affine normalization between channels and arrays"
#
# \description{
#   @get "title".
#
#   This method will both remove curvature in the M vs A plots that are
#   due to an affine transformation of the data. In other words, if there
#   are (small or large) biases in the different (red or green) channels,
#   biases that can be equal too, you will get curvature in the M vs A plots
#   and this type of curvature will be removed by this normalization method.
#
#   Moreover, if you normalize all slides at once (recommended), this 
#   method will also bring the signals on the same scale such that the
#   log-ratios for different slides are comparable. Thus, do not normalize
#   the scale of the log-ratios between slides afterward.
#
#   It is recommended to normalize as many slides as possible in one run.
#   The result is that if creating log-ratios between any channels and any
#   slides, they will contain as little curvature as possible.
#
#   Furthermore, since the relative scale between any two channels on any
#   two slides will be one if one normalizes all slides (and channels) at 
#   once it is possible to add or multiply with the \emph{same} constant 
#   to all channels/arrays without introducing curvature. Thus, it is 
#   easy to rescale the data afterwards as demonstrated in the example.
# }
#
# @synopsis
#
# \arguments{
#  \item{slides}{@vector of slides to be normalized \emph{at once}. 
#   If @NULL, all slides are used.
#   To normalize several slides but slide by slide, call this method
#   once for each slide instead.}
#  \item{groupBy}{@character string or @see "aroma::LayoutGroups" 
#   specifying the groups of spots that are to normalized individually. 
#   If @NULL, all spots are considered.}
#  \item{...}{Additional arguments accepted by 
#   @see "aroma.light::normalizeAffine.matrix".}
# }
#
# \value{
#   Returns a @list containing estimated biases, slopes, and distances
#   to the optimal line for each slide. 
# }
#
# \section{Negative, non-positive, and saturated values}{
#  Affine normalization applies equally well to negative values. Thus,
#  contrary to normalization methods applied to log-ratios, such as curve-fit
#  normalization methods, affine normalization, will not set these to @NA.
#
#  Data points that are saturated in one or more channels are not used 
#  to estimate the normalization function, but they are normalized.
# }
#
# \section{Missing values}{
#  The estimation of the affine normalization function will be made based
#  on only complete non-saturated observations, i.e. observations that
#  contains no @NA values nor saturated values as defined by \code{satSignal}.
# }
#
# \section{Weighted normalization}{
#  Each data point can be assigned a weight in [0,1] specifying how much
#  it should \emph{affect the fitting of the affine normalization function}.
#  Note that here a \emph{data point} is here considered to be the @vector 
#  of values for all channels and all arrays included in the normalization. 
#  For instance, for two-channel data with three arrays, each data point is
#  a vector (R1,G1,R2,G2,R3,G3) of length 6 where R and G represent the
#  red and the green channels.
#
#  Regardless of weights, \emph{all} data points are 
#  \emph{normalized} based on the fitted normalization function.
#
#  Data-point weights are obtained from probe weights, if given.
#  Weights can be set using @seemethod "setProbeWeights".
#  If weights are specified, they will be used.
# }
#
# \section{Within-slide normalization}{
#   This normalization method normalized multiple channels/arrays at once.
#   To normalize array by array, like curve-fit normalization, use a loop,
#   e.g. \code{for (kk in 1:nbrOfSlides(rg)) normalizeAffine(rg, slide=kk)}.
# }
#
# \details{
#  A line is fitted robustly throught the \eqn{(y_R,y_G)} observations
#  using an iterated re-weighted principal component analysis (IWPCA),
#  which minimized the residuals that are orthogonal to the fitted line.
#  Each observation is down-weighted by the inverse of the absolute
#  residuals, i.e. in \eqn{L_1}. 
# }
#
# @author
#
# \references{
#   [1] @include "../incl/BengtssonHossjer_2006.bib.Rdoc" \cr
# }
#
# \examples{
#  @include "../incl/RGData.functions.Rex"
#  @include "../incl/RGData.normalizeAffine.Rex"
# }
#
# \seealso{
#   Internally, the light-weight function 
#   @see "aroma.light::normalizeAffine.matrix" is used.
#
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("normalizeAffine", "RGData", function(this, slides=NULL, groupBy=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Verify the arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'slides'
  slides <- validateArgumentSlides(this, slides=slides);
  
  # Argument: 'groupBy'
  groupBy <- validateArgumentGroupBy(this, groupBy=groupBy);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Prepare the data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Channels to be normalized. Note that to normalize for instance four-
  # channel microarray data, we just have to add the other channels below.
  # /HB 2003-12-29
  channels <- c("R", "G");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Generate data-point weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  weights <- NULL;
  if (hasWeights(this)) {
    hasSignalWeights <- TRUE;
    for (ch in channels)
      hasSignalWeights <- hasSignalWeights && hasSignalWeights(this, channel=ch);

    weights <- matrix(1, nrow=nbrOfSpots(this), ncol=length(slides));

    for (slide in slides) {
      if (hasSignalWeights) {
        w <- matrix(1, nrow=nbrOfSpots(this), ncol=length(channels));
        colnames(w) <- channels;
        for (ch in channels)
          w[,ch] <- getSignalWeights(this, channel=ch, slides=slide);
        # For each array, unite the signal weights into probe weights.
        pw <- MicroarrayWeights$unite(w);
        rm(w);
      } else if (hasProbeWeights(this)) {
        # The data-point weight is the square root of the average squares.
        pw <- getProbeWeights(this, slides=slide);
      }

      weights[,slide] <- pw;
      rm(pw);
    }

    # Unite all probe weights into data points weights for all arrays.
    weights <- MicroarrayWeights$unite(weights);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Fit the model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fits <- list();
  # For all groups of spots...
  for (spots in groupBy) {
    # Put all slides and arrays in the same matrix such the data from the
    # first array comes first, then the second array and so on, e.g.
    # c(R1,G1,R2,G2,R3,G3) for 3 two-color arrays.
    X <- NULL;
    for (slide in slides) {
      for (ch in channels) {
        X <- cbind(X, this[[ch]][spots,slide]);
      }
    }
    colnames(X) <- rep(channels, times=length(slides));

    # Normalize the data
    X <- normalizeAffine(X, weights=weights[spots], ...);
    
    fit <- attr(X, "modelFit");
#    fit <- AffineModelFit(t0=fit$t0, t=fit$t, a=fit$a, b=fit$b, adiag=fit$adiag, eigen=fit$U, y=fit$X, converged=fit$converged, nbrOfIterations=fit$niter); 
    fits <- c(fits, list(fit));
    mostattributes(X) <- NULL;
    
    # Store the normalized data
    for (ch in channels) {
      cols <- (colnames(X) == ch);
      this[[ch]][spots,slides] <- X[,cols];
    }

    # Clean up
    rm(fit, X);
  } # for (spots in groupBy) 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Post-processing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Data object has changed, clear cached values.
  clearCache(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return the parameter estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  invisible(fits);
}) # normalizeAffine()





#########################################################################/**
# @RdocMethod normalizeCurveFit
# @alias normalizeLowess
# @alias normalizeLoess
# @alias normalizeSpline
# @alias normalizeRobustSpline
#
# @title "Within-slide intensity-dependent normalization in (A,M)"
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
# \value{
#   Returns itself invisibly.
# }
#
# \description{
#   @get "title". It normalizes pairs of channels by estimating a smooth 
#   log-ratio intensity-dependent curve.
#
#   Because of they are so well known by their de facto names, the methods
#   \code{normalizeLowess(...)} and \code{normalizeLoess(...)} are aliases 
#   for \code{normalizeCurveFit(..., method="lowess")} and
#   \code{normalizeCurveFit(..., method="loess")}, respectively.
# }
#
# \section{Negative, non-positive, and saturated values}{
#  Non-positive values are set to not-a-number (@NaN).
#  Data points that are saturated in one or more channels are not used 
#  to estimate the normalization function, but they are normalized.
# }
#
# \section{Missing values}{
#  The estimation of the affine normalization function will be made based
#  on only complete non-saturated observations, i.e. observations that
#  contains no @NA values nor saturated values as defined by \code{satSignal}.
# }
#
# \section{Weighted normalization}{
#  Each data point can be assigned a weight in [0,1] specifying how much
#  it should \emph{affect the fitting of the curve-fit normalization function}.
#  Note that here a \emph{data point} is here considered to be the pair
#  of values in the two channels to be normalized. 
#  For instance, for two-channel data, a data point is the pair (R,G).
#
#  Regardless of weights, \emph{all} data points are \emph{normalized} based 
#  on the fitted normalization function.
#
#  Weights can be set using @seemethod "setProbeWeights".
#  If weights are specified, they will be used.
# }
#
# \examples{
#  @include "../incl/RGData.functions.Rex"
#  @include "../incl/RGData.normalizeCurveFit.Rex"
# }
#
# @author
#
# \seealso{
#   Internally, the light-weight function 
#   @see "aroma.light::normalizeCurveFit.matrix" is used.
#
#   @seeclass
# }
#*/#########################################################################
setMethodS3("normalizeCurveFit", "RGData", function(this, slides=NULL, groupBy=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'slides'
  slides <- validateArgumentSlides(this, slides=slides);

  # Argument: 'groupBy'
  spots <- validateArgumentGroupBy(this, groupBy=groupBy);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (slide in slides) {
    Xall <- cbind(this$R[,slide], this$G[,slide]);
    wall <- getProbeWeights(this, slides=slide);
    for (kk in seq(along=spots)) {
      spot <- spots[[kk]];
      X <- Xall[spot,];
      w <- wall[spot,];
      X <- normalizeCurveFit(X, weights=w, ...);
      Xall[spot,] <- X;
    } # for (spot ...)
    rm(X,w);
    this$R[,slide] <- Xall[,1];
    this$G[,slide] <- Xall[,2];
  } # for (slide ...)
  rm(Xall,wall);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Post-processing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Data object has changed, clear cached values.
  clearCache(this);

  invisible(this);
}) # normalizeCurveFit();


setMethodS3("normalizeLowess", "RGData", function(this, ...) {
  normalizeCurveFit(this, method="lowess", ...);
})

setMethodS3("normalizeLoess", "RGData", function(this, ...) {
  normalizeCurveFit(this, method="loess", ...);
})

setMethodS3("normalizeSpline", "RGData", function(this, ...) {
  normalizeCurveFit(this, method="spline", ...);
})

setMethodS3("normalizeRobustSpline", "RGData", function(this, ...) {
  normalizeCurveFit(this, method="robust.spline", ...);
})






#########################################################################/**
# @RdocMethod calibrateMultiscan
#
# \encoding{latin1}
# 
# @title "Calibrates multiple re-scanned images based on an affine model"
#
# \description{
#   @get "title".
#
#   Each channel is calibrated seperately.
# }
#
# @synopsis
#
# \arguments{
#  \item{slides}{Slides to be used in the fit \emph{and} that are to be 
#   calibrated. If @NULL, all slides are considered.}
#   \item{channels}{@character string specifying which channels to be
#    calibrated. If @NULL, all channels are calibrated.}
#   \item{groupBy}{@character string or @see "aroma::LayoutGroups" 
#     specifying the groups of spots that are to calibrated individually. 
#     If @NULL, all spots are considered.}
#  \item{...}{Additional arguments accepted by
#     @see "aroma.light::calibrateMultiscan.matrix". 
#     Its help page is very USEFUL.}
# }
#
# \value{
#   Returns a @list containing one element for each calibrated channel.
#   Each channel element contains parameter estimates either directly
#   as global estimates or as a @list consisting estimates for each
#   group (as defined \code{groupBy}).
# }
#
# \section{Weighted calibration}{
#  Each data point can be assigned a weight in [0,1] specifying how much
#  it should \emph{affect the fitting of the calibration function}.
#  Note that here a \emph{data point} is here considered to be the @vector 
#  of values from all scans ("slides").
#
#  Regardless of weights, \emph{all} data points are \emph{calibrated} based 
#  on the fitted calibrated function.
#
#  Data-point weights are obtained from probe weights, if given.
#  Weights can be set using @seemethod "setProbeWeights".
#  If weights are specified, they will be used.
#  Currently it is not possible to set different in different channels.
# }
#
# \details{
#  Fitting is done by iterated re-weighted principal component analysis
#  (IWPCA).
# }
#
# @author
#
# \references{
#   [1] @include "../incl/BengtssonH_etal_2004.bib.Rdoc" \cr
# }
#
# \examples{
#  @include "../incl/RGData.functions.Rex"
#  @include "../incl/RGData.calibrateMultiscan.Rex"
# }
#
# \seealso{
#   @see "aroma.light::calibrateMultiscan.matrix".
#   @see "normalizeAffine.RGData".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("calibrateMultiscan", "RGData", function(this, slides=NULL, channels=NULL, groupBy=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'slides'
  slides <- validateArgumentSlides(this, slides=slides);
  
  # Argument: 'channels'
  channels <- c("R", "G");
  if (is.null(channels) || identical(channels, "all")) {
    channels <- c("R", "G");
  } else {
    channels <- as.character(channels);
  }

  # Argument: 'groupBy' 
  groupBy <- validateArgumentGroupBy(this, groupBy=groupBy);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Generate data-point weights from 
  # i) signal weights, or ii) probe weights.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  weights <- NULL;
  if (hasWeights(this)) {
    weights <- list();
    for (ch in channels) {
      if (hasSignalWeights(this, channel=ch)) {
        w <- getSignalWeights(this, channel=ch, slides=slides);
        # For each channel, unite the signal weights from several scans 
        # into data point weights.
        w <- MicroarrayWeights$unite(w);
      } else if (hasProbeWeights(this)) {
        # The data-point weight is the square root of the average squares.
        w <- getProbeWeights(this, slides=slides);
        w <- MicroarrayWeights$unite(w);
      }
      weights[[ch]] <- w;
      rm(w);
    }
  }
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Fit the model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fits <- list();
  
  # Calibrate each channel seperately
  for (channel in channels) {
    fits[[channel]] <- list();
    
    # Calibrate each spot group
    for (spots in groupBy) {
      X <- this[[channel]][spots, slides];
      w <- weights[[channel]][spots];

      # Multiscan calibration
      nbrOfSlidesBefore <- ncol(X);
      X <- calibrateMultiscan(X, weights=w, ...);
      nbrOfSlidesAfter <- ncol(X);

      fit <- attr(X, "modelFit");
      fits[[channel]] <- c(fits[[channel]], fit);
      mostattributes(X) <- NULL;
      
      this[[channel]][spots,slides] <- X;
    } # for (channel ...)
  } # for (spots in groupBy)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Post-processing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!is.null(weights)) {
    # Set new signal weights
    for (ch in channels)
      setSignalWeights(this, channel=ch, weights=weights[[ch]], slides=slides);
  }

  # Exclude slides
  if (nbrOfSlidesAfter < nbrOfSlidesBefore) {
    excl <- slides[seq(from=nbrOfSlidesAfter+1, to=nbrOfSlidesBefore)];
    removeSlides(this, slides=excl);
  }

  # Data object has changed, clear cached values.
  clearCache(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Return the parameter estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  invisible(fits);
}) # calibrateMultiscan()




setMethodS3("getCalibratedMultiscan", "RGData", function(this, ...) {
  rgC <- clone(this);
  calibrateMultiscan(rgC, project=TRUE, ...);
  # After projecting the data onto the theoretical line in the model all
  # 'slides' are identical. Return the first one...
  as.RGData(rgC, slide=1);
}) # getCalibratedMultiscan()




############################################################################
# HISTORY:
# 2005-02-12
# o Updated Rdoc comments to link to calibrateMultiscan(), 
#   normalizeAffine(), and normalizeCurveFit() in aroma.light.
# 2005-02-08
# o Now calibrateMultiscan() and normalizeAffine() supports signal weights. 
#   This is useful since we can downweight spots at high intensities, 
#   because they might be (partly) saturated. Moreover, calibrateMultiscan()
#   sets signal weights based on united previous signal weights.
# 2005-02-04
# o Remove getAverageMultiscan() and averageMultiscan(). Now (again) taken
#   care of by calibrateMultiscan().
# 2005-02-02
# o Changed the order how normalizeAffine() puts the channels. It is easier
#   to return estimates if it is done (R1,G1,R2,G2,...).
# o Excluded fitIWPCA() and fitMultiscanAffine() from RGData() since they
#   are now in the matrix class.
# o Now weights in normalizeAffine() are calculated from probe weights.
# 2005-02-01
# o Added normalizeCurveFit(), which now utilizing the method with the
#   same name, but for the matrix class.
# 2004-05-09
# o calibrateMultiscan() is now using the calibrateMultiscan.matrix()
#   method. Similar for the normalizeAffine() method.
# o BUG FIX: The background transform in calibrateMultiscan() did not
#   transform the matrix before subtracting. This made no difference when
#   using the default constraint=="diagonal", but would introduce strange
#   biases in others would be used.
# 2004-04-26
# o Had to add argument 'w' here so that the weights are also "resampled"
#   together with the data points in the internal bootstrap.
# o Reverted back to estimate b in fitIWPCA() such that b[1] == 1.
#   Otherwise estimates such as printtip and spatial are really odd;
#   for some groups b[1] was equal to one and sometimes b[2] was one.
# 2004-04-18
# o Added the ... argument that is passed down all the way to the iwpca()
#   method. This allows for instance a weight argument to normalizeAffine()
#   and calibrateMultiscan(), and to the internal fitMultiscanAffine() and
#   fitIWPCA().
# 2004-03-10
# o Updated fitIWPCA() to make use of iwpca(), which in turn makes use of
#   our own wpca() instead of as before acp() in multidim.
# 2004-03-01
# o Added argument 'channels' to calibrateMultiscan().
# o Made argument 'project' of  calibrateMultiscan FALSE by default.
# 2004-02-17
# o Moved from the aroma.affine package to the aroma package.
# 2004-01-12
# o Added argument 'groupBy' to normalizeAffine() and calibrateMultiscan(), 
#   which means that they now both support print-tip, plate etc
#   normalization too.
# 2004-01-07
# o BUG FIX: Several occurances of foo <- match.arg("foo") - removed '"'s.
# o Added the constraints "diagonal" and "max" to fitIWPCA().
# o Added the argument 'constraint' to normalizeAffine().
# 2004-01-02
# o Made fitIWPCA() call gc() before every estimation.
# 2004-01-01
# o Made fitIWPCA() work with boot().
# 2003-12-30
# o Added removeSlides() and keepSlides()
# 2003-12-29
# o Making use of distanceBetweenLines() in R.basic now. Thus, the code is
#   again made much cleaner. Thanks to this, we also do not have to treat
#   the two-dimensional situation seperately.
# 2003-12-28
# o Thanks to iwpca() we can put the code for estimating the parameters
#   whe K == 2 inside function too.
# o Extracted the code for robust PCA into iwpca() now in the R.basic 
#   package. This makes this code much cleaner.
# o fitIWPCA() replaces former fitMultiIWPCA() and fitPairIWPCA().
#   Now the method returns an object of class AffineMultiscanModelFit.
# o fitMultiscanAffine(): Adjusted to make use of the new
#   AffineMultiscanModelFit class. Renamed parameter 'k' to 'b' as in
#   the paper.
# o calibrateMultiscan(): Updated code to make use of the new
#   AffineMultiscanModelFit class.
# o Updated normalizeAffine() to make use of RGData$fitIWPCA().
# 2003-12-27
# o Added Rdoc comments and much more source code comments.
# 2003-12-10
# o Now the returned parameter vector 'k' is normalized such that k[1] == 1.
#   This is also is line with fitPairIWPCA().
# o fitMultiscanAffine(): Added first support for bootstrapping the 
#   parameters and also the 'control' argument.
# 2003-11-09
# o Added getWithinChannelPairs().
# 2003-12-06
# o Added the argument 'project' to calibrateMultiscan() to make it explicit
#   that if we project the data onto the fitted line, all scans will get the
#   identical signals afterwards (simply because they are projected onto the
#   same point!).
# o BUG FIX: Of course, normalizing between slides, which is based on the
#   *assumption* that most genes are *approximately* non-differentially
#   expressed, does *not* justify that we can project the data (in 2D) onto
#   the fitted line as we do for single-channel calibrate multiscan. This
#   will of course make the signal in both channels identical, which should
#   only be true for truely non-differentially expressed genes. Thus, no
#   projection algorithm for normalizeAffine(). In calibrateMultiscan() we
#   *know* that the signals should be identical (or at least almost).
# 2003-11-10
# o Added getCalibratedMultiscan().
# o Now the method calibrateMultiscan() projects the data onto the translated
#   line L' that goes through a*(1,1,...1) to obtain ytilde. These points are
#   the used to get xtilde = (ytilde-a)/b which are our estimates.
# 2003-11-07
# o The "translation error" in a, and the eigenvectors and the center of
#   the principal components are now also returned. Based on Ola Hössjer's
#   suggestion to be used in our factor analysis model.
# 2003-11-05
# o Added fitMultiscanSpatial().
# 2003-10-29
# o Added the arguments 'mad=TRUE' and 'sd=FALSE' to averageMultiscan().
# 2003-10-26
# o Extracted the "fit" part of calibrateMultiscan() to fitMultiscanAffine()
#   to be able to test different fits. This is only for internal use.
# o Extract the fitting of the affine model to its own method, which
#   calibrateMultiscan() is calling and then just calibrating.
# o Added argument 'slides' to calibrateMultiscan() too. This will allow
#   too load all slides at once an calibrate each slide seperately. It also
#   allows you to do pairwise calibrations to for instance compare different
#   scan pairs for for instance bleach effects etc.
# 2003-10-25
# o Renamed to calibrateMultiscan() according to my talk.
# 2003-10-19
# o Extracted from normalizeAffine.R as these functions are used by both
#   normalizeAffine() and normalizeAffineMultiscan().
# o Created normalizeAffine() from former normalizeAffineMultiscan().
# 2003-10-18
# o Created calibrateMultiscan() from previous PMT analysis scripts.
############################################################################
