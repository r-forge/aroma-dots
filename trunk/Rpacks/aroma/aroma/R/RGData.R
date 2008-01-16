#########################################################################/**
# @RdocClass RGData
#
# @title "The RGData class"
#
# \description{
#  @classhierarchy
#
#  Creates a new \code{RGData} object. 
#  The philosophy behind this data structure is to think about the data in the form of the signals in one channel (R) versus the signals in the other channel (G). 
#  This is in contrast to the idea of the @see "MAData" structure, which thinks about the data as the log ratios (M) and log intensites (A) for the spot signals. 
# }
#
# @synopsis
#
# \arguments{
#   \item{R,G}{A NxM @matrix containing (non-logged) signals of the red
#    (green) channel, where N is the number of spots on each slide and M
#    is the number of slides in this data set.}
#   \item{layout}{A @see "Layout" object specifying the spot layout of the
#    slides in this data set.}
#   \item{extras}{Private argument. Do not use.}
# }
#
# \section{Fields and Methods}{
#  \bold{Fields}
#  \tabular{rll}{
#    \tab \code{R} \tab The signal for channel R (non-logged). \cr
#    \tab \code{G} \tab The signal for channel G (non-logged). \cr
#  }
#
#  @allmethods "public"
# }
#
# \details{
#   The mapping between M and A, and R and G is a one-to-one function.
#   Given the signal R and G for the R and the G channels you get the
#   M and the A values by:
#   \deqn{
#     M = \log_2\frac{R}{G},\quad 
#     A = \log_2\sqrt{R{\cdot}G} = \frac{1}{2}\log_2 R{\cdot}G,
#   }{
#     M = log2(R/G), A = log2(sqrt(R*G)) = 1/2*log2(R*G),
#   }
#   which in [R] can be done by \code{ma <- as.MAData(rg)}. The reverse function, i.e. going back to the R and the G is:
#   \deqn{
#     R = \sqrt{2^{2A+M}},\quad G = \sqrt{2^{2A-M}}
#   }{
#     R = sqrt(2^(2A+M)), G = sqrt(2^(2A-M))
#   }
#   which in [R] can be done by \code{rg <- as.RGData(rg)}.
#
#   Note that if the signal in one or both channels is non-positive, 
#   the log-transform will make these values undefined, that is, set
#   them to @NA. When going back to (G,R) from (A,M) these values
#   will remain @NA. 
# }
#
# \note{
#  There are several functions that returns an object of this class, and
#  it is only in very special cases that you actually have to create one
#  yourself.
# }
#
# @author
#
# \examples{
#   # Create a raw data object from the preexisting example data in
#   # the sma package.
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#
#   # Get the signal (here by default non-background corrected)
#   ma <- getSignal(raw)
#
#   # Transform (M,A) into (R,G)
#   rg <- as.RGData(ma)
# }
#*/#########################################################################
setConstructorS3("RGData", function(R=NULL, G=NULL, layout=NULL, extras=list()) {
  if (!is.null(R)) R <- SpotSlideArray(R, layout=layout);
  if (!is.null(G)) G <- SpotSlideArray(G, layout=layout);
  this <- extend(MicroarrayData(layout=layout, extras=extras), "RGData",
    R = R,
    G = G
  )
  this$.fieldNames <- c("R", "G");
  
  # Sets the default labels:
  setLabel(this, "R", expression(R));
  setLabel(this, "G", expression(G));

  this;
})


setMethodS3("as.character", "RGData", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- "";
  s <- paste(sep="",s,"R ",   com.braju.sma.dimStr(this$R));
  s <- paste(sep="",s,", G ", com.braju.sma.dimStr(this$G));
  if (hasLayout(this))
    s <- paste(sep="", s, " ", getLayout(this));
  if (hasWeights(this)) {
    w <- sapply(this$weights, FUN=is.null);
    s <- paste(s, " Weights (", paste(names(w), collapse=", "), 
                  ") are specified.", sep="");
  }
  if (!is.null(this$.cache) && length(this$.cache) > 0)
    s <- paste(sep="",s,", cached values for ",
               paste(names(this$.cache), collapse=", "), " exist");
  s;
})


setMethodS3("getChannelNames", "RGData", function(this) {
  c("R", "G");
})


#########################################################################/**
# @RdocMethod swapDyes
#
# @title "Swap dyes of one or many slides"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \value{
#   Return itself.
# }
#
# \arguments{
#   \item{slides}{A @vector of slides to be dye swapped. If @NULL, all
#   slides are considered.}
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#
#   ma <- getSignal(raw)
#   rg <- as.RGData(ma)
#
#   # Dye swap every other slide.
#   swapDyes(rg, slides=c(4,5,6))
#
#   layout(matrix(1:6, nrow=2, ncol=3, byrow=TRUE));
#   for (k in 1:6)
#     plot(rg, slide=k)
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("swapDyes", "RGData", function(this, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  
  tmp <- this$R[,slides];
  this$R[,slides] <- this$G[,slides];  
  this$G[,slides] <- tmp;

  clearCache(this); 
  
  invisible(this);
})




############################################################################
############################################################################
## 
##  DATA STRUCTURAL METHODS
## 
############################################################################
############################################################################

setMethodS3("as.RGData", "ANY", function(obj, ...) {
  RGData(obj, ...)
})

setMethodS3("as.RGData", "RawData", function(this, ...) {
  as.RGData.RGData(this, ...);
})


setMethodS3("as.RGData", "RGData", function(this, slides=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);

  R <- this$R[,slides];
  G <- this$G[,slides];
  layout <- getLayout(this);
  extras <- this$.extras;
  res <- RGData(R=R, G=G, layout=layout, extras=extras);
  res$weights <- getWeights(this, slides=slides);

  res;
});




#########################################################################/**
# @RdocMethod getLogRatios
# @alias "getM"
#
# @title "Calculates the log-ratios (M values)"
#
# \description{
#   @get "title" for each specified array. 
#   Log-base 2 (two) is used.
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{A @vector of @integers indicating which slides to be
#     considered. If @NULL, all slides are considered.}
# }
#
# \value{
#   Returns a NxK @matrix where N is the number of spots and K is the 
#   number of (specified) arrays.
# }
#
# \seealso{
#   @seemethod "getLogIntensities".
#   @seeclass
# }
#
# @author
#*/#########################################################################
setMethodS3("getLogRatios", "RGData", function(this, slides=NULL) {
  getM(this, slides=slides);
})

setMethodS3("getM", "RGData", function(this, slides=NULL) {
  log.na <- function(x, ...) {
    ifelse(x > 0, log(x, ...), NA)
  }

  slides <- validateArgumentSlides(this, slides=slides);

  R <- this$R[,slides];
  G <- this$G[,slides];

  # To prevent R*G to produce "integer overflow" make sure it is calculated
  # in double precision!
  R <- matrix(as.double(R), ncol=length(slides));

  # Preserve the slide names
  colnames(R) <- getSlideNames(this, slides=slides);

  # If any of the two channel contains negative signals, make them NA's
  R[R < 0] <- NA;
  G[G < 0] <- NA;

  SpotSlideArray(log.na(R/G, 2), layout=getLayout(this));
}, protected=TRUE);




#########################################################################/**
# @RdocMethod getLogIntensities
# @alias "getA"
#
# @title "Calculates the log-intensitites (A values)"
#
# \description{
#   @get "title" for each specified array.
#   Log-base 2 (two) is used.
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{A @vector of @integers indicating which slides to be
#     considered. If @NULL, all slides are considered.}
# }
#
# \value{
#   Returns a NxK @matrix where N is the number of spots and K is the 
#   number of (specified) arrays.
# }
#
# \seealso{
#   @seemethod "getLogRatios".
#   @seeclass
# }
#
# @author
#*/#########################################################################
setMethodS3("getLogIntensities", "RGData", function(this, slides=NULL) {
  getA(this, slides=slides);
})

setMethodS3("getA", "RGData", function(this, slides=NULL) {
  log.na <- function(x, ...) {
    ifelse(x > 0, log(x, ...), NA)
  }

  slides <- validateArgumentSlides(this, slides=slides);

  R <- this$R[,slides];
  G <- this$G[,slides];

  # To prevent R*G to produce "integer overflow" make sure it is calculated
  # in double precision!
  R <- matrix(as.double(R), ncol=length(slides));

  # Preserve the slide names
  colnames(R) <- getSlideNames(this, slides=slides);

  # If any of the two channel contains negative signals, make them NA's
  R[R < 0] <- NA;
  G[G < 0] <- NA;

  SpotSlideArray(log.na(R*G, 2)/2, layout=getLayout(this));
}, protected=TRUE);




#########################################################################/**
# @RdocMethod as.MAData
#
# @title "Transform from the red and green intensities into log ratios between them and the log product of them"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{A @vector of @integers indicating which slides to transformed. If @NULL, all slides are transformed.}
# }
#
# \value{
#   Returns a @see "MAData" object.
# }
#
# \examples{
#   SMA$loadData(c("mouse.data", "mouse.setup"))
#   raw <- RawData(mouse.data, layout=as.Layout(mouse.setup))
#   rg <- getSignal(raw)
#   ma <- as.MAData(rg)
# }
#
# \seealso{
#   @see "MAData.as.RGData".
#   @seeclass
# }
#
# @author
#*/#########################################################################
setMethodS3("as.MAData", "RGData", function(this, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  R <- this$R[,slides];
  G <- this$G[,slides];

  # To prevent R*G to produce "integer overflow" make sure it is calculated
  # in double precision!
  R <- matrix(as.double(R), ncol=length(slides));

  # Preserve the slide names
  colnames(R) <- getSlideNames(this, slides=slides);

  # If any of the two channel contains negative signals, make them NA's
  R[R < 0] <- NA;
  G[G < 0] <- NA;
  MAData(M=log.na(R/G, 2), A=log.na(R*G, 2)/2, layout=getLayout(this), 
                                                extras=this$.extras)
})






setMethodS3("as.RawData", "RGData", function(this, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  R <- this$R[,slides];
  G <- this$G[,slides];
  Rb <- matrix(0, ncol=length(slides));
  Gb <- Rb;

  RawData(R=R, G=G, Rb=Rb, Gb=Gb, layout=getLayout(this), 
                                                extras=this$.extras)
})
            


############################################################################
# @RdocMethod getWithinChannelPairs
#
# @title "Gets all possible slide combinations of one of the channels"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{channel}{The channel whose signals should be returned.}
#   \item{slides}{Subset of slides to be used. If @NULL, all slides
#      are considered.}
# }
#
# \value{
#   Returns the background signal as a @see "RGData" object.
# }
#
# \examples{\dontrun{See help(RawData.getWithinChannelPairs) for an example.}}
#
# @author
# 
# \seealso{
#   See also \code{getSlidePairs()} in the @see "MicroarrayData@ class,
#   which is used internally.
#   @seeclass
# }
############################################################################
setMethodS3("getWithinChannelPairs", "RGData", function(this, channel, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  pairs <- getSlidePairs(this, slides=slides);

  X1 <- this[[channel]][,pairs[1,]];
  X2 <- this[[channel]][,pairs[2,]];

  RGData(R=X1,G=X2, layout=getLayout(this));
})



############################################################################
############################################################################
## 
##  PLOTTING & GRAPHICAL METHODS
## 
############################################################################
############################################################################


#########################################################################/**
# @RdocMethod getColors
# 
# @title "Generates red to green colors for each of the specified spots"
#
# \description{
#   @get "title", which can
#  be passed to the \code{col} argument in most plot functions.
# }
#
# @synopsis
#
# \arguments{
#  \item{palette}{The color palette to be used.}
#  \item{M.range}{The range of the M=log2(R/G). Values outside this range will be saturated to the extreme green and red, respectively.}
#  \item{A.range}{The range of the A=1/2*log2(R*G). Values outside this range will be saturated to the dark and light, respectively.}
# }
#
# \value{Returns nothing.}
#
# @author
#
# \examples{
#   SMA$loadData(c("mouse.data", "mouse.setup"))
#   raw <- RawData(mouse.data, layout=as.Layout(mouse.setup))
#   ma <- getSignal(raw)
#   rg <- as.RGData(ma)
#   col <- getColors(rg)
# }
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getColors", "RGData", function(this, what=c("A","M"), slide=1, include=NULL, exclude=NULL, palette=NULL, range=matrix(c(0,16, -2,2), nrow=2), ...) {
  if (length(what) == 2) {
    ma <- as.MAData(this);
    what <- c("A", "M");
    ma$getColors(what=what, slide=slide, include=include, exclude=exclude, 
                                        palette=palette, range=range);
  } else {
    what <- what[1];
    x <- this[[what]];
    x <- log(x, 2);
    x[is.na(x)] <- 0;

    if (is.null(palette)) palette <- "auto";
    palettes <- c("grayscale", "auto", "redgreen", "redscale", "greenscale");
    idx <- match(palette, palettes);
    if (is.null(idx))
      throw("Unknown palette (\"", palette, "\").");

    if (idx == 1) {
      Colors$getGray(x, x.range=range);
    } else {
      channel <- substring(what,1,1);
      dimensions <- c("", channel, channel, "R", "G");
      dim <- dimensions[idx];
      Colors$getRGB(x, x.range=range, dim=dim);
    }
  }
});
            


setMethodS3("plotXY", "RGData", function(this, what=c("G", "R"), xlog=2, ylog=xlog, xlim=NULL, ylim=xlim, ...) {
  plotXY.MicroarrayData(this, what=what, xlog=xlog, ylog=ylog, xlim=xlim, ylim=ylim, ...);
})
            

setMethodS3("plotSpatial", "RGData", function(this, what="R", col="auto", ...) {
  if (!is.null(col)) {
    if (col == "auto") {
      if (what == "R") col <- "redscale"
      else if (what == "G") col <- "greenscale"
      else col <- NULL;
    } else {
      col <- "grayscale";
    }
  }
  plotSpatial.MicroarrayData(this, what=what, col=col, ...);
})
            

setMethodS3("boxplot", "RGData", function(x, what="R", ...) {
  # To please R CMD check...
  this <- x;

  boxplot.MicroarrayData(this, what=what, ...);
})
            

setMethodS3("plot", "RGData", function(x, what="RvsG", ...) {
  # To please R CMD check...
  this <- x;

  plot.MicroarrayData(this, what=what, ...);
})
            



############################################################################
############################################################################
## 
##  STATISTICAL METHODS
## 
############################################################################
############################################################################

setMethodS3("range", "RGData", function(this, what="M", slide=NULL, na.rm=TRUE, inf.rm=TRUE, ...) {
  if (is.element(what, names(this))) {
    X <- this[[what]];
    if (!is.null(slide)) {
      if (!is.matrix(X))
        throw("The field that 'what' refers to is not a applicable field.");
      X <- X[,slide];
    }
    if (inf.rm) X[is.infinite(X)] <- NA;
    range(X, na.rm=na.rm);
  } else
    range.MicroarrayData(this, what=what, na.rm=na.rm, inf.rm=inf.rm, ...);
})



setMethodS3("log", "RGData", function(x, base=exp(1), ...) {
  # To please R CMD check...
  this <- x;

  res <- clone(this);
  res$R <- log(this$R, base=base);
  res$G <- log(this$G, base=base);
  res;
}, protected=TRUE);



setMethodS3("power", "RGData", function(this, base=exp(1)) {
  res <- clone(this);
  if (base == exp(1)) {
    res$R <- exp(this$R);    # More efficient?
    res$G <- exp(this$G);
  } else {
    res$R <- base^(this$R);
    res$G <- base^(this$G);
  }
  res;
}, protected=TRUE);

            
#########################################################################/**
# @RdocMethod mean
#
# @title "Genewise Average Mean for channel R and G"
#
# \description{
#   Computes the genewise average mean across all slides in the
#   @see "RGData" object.
# }
#
# @synopsis
#
# \arguments{
#  \item{inf.rm}{A @logical value indicating whether @Inf values should be
#         stripped before the computation proceeds.}
#  \item{na.rm}{A @logical value indicating whether @NA values should be
#         stripped before the computation proceeds.}
# }
#
# \value{Returns a new @see "RGData" object containing the result.}
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("mean", "RGData", function(x, inf.rm=TRUE, na.rm=TRUE, ...) {
  # To please R CMD check...
  this <- x;

  res <- clone(this);
  for (name in c("R", "G")) {
    data <- this[[name]];
    ncol <- ncol(data);
    nrow <- nrow(data);
    if (inf.rm)
      data <- matrix(data[!is.infinite(data)], ncol=ncol, nrow=nrow);
    res[[name]] <- as.matrix(apply(data, MARGIN=1, mean, na.rm=na.rm, ...));
  }
  res;
});
            



#########################################################################/**
# @RdocMethod var
#
# @title "Genewise Variance for channel R and G"
#
# \description{
#   Computes the genewise variance across all slides in the
#   @see "RGData" object.
# }
#
# @synopsis
#
# \arguments{
#  \item{inf.rm}{A @logical value indicating whether @Inf values should be
#         stripped before the computation proceeds.}
#  \item{na.rm}{A @logical value indicating whether @NA values should be
#         stripped before the computation proceeds.}
# }
#
# \value{Returns a new @see "RGData" object containing the result.}
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("var", "RGData", function(this, inf.rm=TRUE, na.rm=TRUE, ...) {
  res <- clone(this);

  for (name in c("R", "G")) {
    data <- this[[name]];
    ncol <- ncol(data);
    nrow <- nrow(data);
    if (inf.rm)
      data <- matrix(data[!is.infinite(data)], ncol=ncol, nrow=nrow);
    res[[name]] <- as.matrix(apply(data, MARGIN=1, var, na.rm=na.rm, ...));
  }
  res;
})
            



setMethodS3("read", "RGData", function(this, filename, path=NULL, layout=NULL, verbose=FALSE) {
  fields <- c("R", "G");
  res <- MicroarrayData$readToList(filename, path=path,
                                   reqFields=fields, verbose=verbose);
  
  # Create the MAData object.
  RGData(R=res$R, G=res$G, layout=layout)
}, static=TRUE, trial=TRUE);



setMethodS3("normalizeGenewise", "RGData", function(this, fields=c("R", "G"), bias=c(10,10), scale=1, ...) {
  NextMethod("normalizeGenewise", this, fields=fields, bias=bias, scale=scale, ...);
})


setMethodS3("normalizeQuantile", "RGData", function(this, 
                                                  fields=c("R", "G"), ...) {
  normalizeQuantile(this, fields=fields, ...);
}) # normalizeQuantile()





#########################################################################/**
# @RdocMethod shift
#
# @title "Shift the log-ratios, log-intensities or the raw signal"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{M,A,R,G}{A @numeric or @function specifying the shift to be 
#    applied to the log-ratios, the log-intensities, the red signals, 
#    and/or the green signals.
#    If more than one of these are shifted at the same time, they are
#    shifted in the order \code{M}, \code{A}, \code{R} and \code{G}.
#    A @numeric specify the amount of shift. 
#    If a @function, e.g. \code{min()}, is used, then the amount of shift
#    is the value returned by that function when all \emph{finite} values
#    are passed to that function, e.g. \code{min(x[is.finite(x)])}. In 
#    other words, @NA's etc are automatically taken care of.
#   }
#   \item{slides}{Slides to be shifted. If @NULL, all slides are shifted.}
# }
#
# \value{
#   Returns nothing.
# }
#
# \examples{\dontrun{See shift() for the MAData class for an example.}}
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("shift", "RGData", function(this, M=NULL, A=NULL, R=NULL, G=NULL, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  foo <- function(x, X) {
    if (is.function(X)) {
      x <- x - X(x[is.finite(x)]);
    } else {
      x <- x + X;
    }
  }

  for (slide in slides) {
    # i) shift R and G
    if (!is.null(R))
      this$R[,slide] <- foo(this$R[,slide], R);
    if (!is.null(G))
      this$G[,slide] <- foo(this$G[,slide], G);

    # ii) shift M and A
    if (!is.null(M) || !is.null(A)) {
      r <- this$R[,slide];
      g <- this$G[,slide];

      m <- log(r/g, base=2);
      a <- 1/2*log(r*g, base=2);

      if (!is.null(M))
  	m <- foo(m, M);
      if (!is.null(A))
  	a <- foo(a, A);

      r <- sqrt(2^(2 * a + m));
      g <- sqrt(2^(2 * a - m));
      rm(m,a);
      this$R[,slide] <- r;
      this$G[,slide] <- g;
    }
  } # for (slide in slides)

  clearCache(this); 
})



setMethodS3("findChannelBiasDifferences", "RGData", function(this, slides=NULL, method=c("robust-pca", "pca", "loess-extrapolate", "loess-lm"), maxIter=30, visualize=FALSE, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  method <- match.arg(method);

  res <- list();

  for (slide in slides) {
    R <- this$R[,slide];
    G <- this$G[,slide];

    if (visualize)
      plot(G, R, pch=176, xlim=c(0,2e3), ylim=c(0,2e3));

    if (method == "pca") {
      if (!require("mva") && !require("stats")) {
        throw("Failed to locate princomp(). Package not loaded: mva or stats");
      }
      fit <- princomp(~R+G, scores=TRUE);
      A <- unclass(loadings(fit));
      x0 <- fit$center;
      mR <- x0[2]-A[2,1]/A[1,1]*x0[1];
      kRG <- A[2,1]/A[1,1];
      offsetR <- c(G0=0, R0=mR);
      offsetG <- c(G0=mR, R0=0);
      offset <- rbind(offsetR, offsetG);
    } else if (method == "robust-pca") {
      acp <- NULL; rm(acp);  # Dummy to please R CMD check.
      require("multidim") || throw("Package 'multidim' not found.");
      reps <- 0.01;
      xy <- cbind(R,G);
      w <- rep(1, length=nrow(xy));
      for (iter in 1:maxIter) {
        fit <- acp(xy, wt=w, reduc=FALSE);
        w <- 1/abs(fit$cmpr[,2]+reps);
      }
      A <- fit$vectors;
      x0 <- fit$moy;
      mR <- x0[2]-A[2,1]/A[1,1]*x0[1];
      kRG <- A[2,1]/A[1,1];
      offsetR <- c(G0=0, R0=mR);
      offsetG <- c(G0=mR, R0=0);
      offset <- rbind(offsetR, offsetG);
    } else if (method == "loess-extrapolate") {
      #  i) Fit a loess curve through the data.
      fit <- loess(R ~ G, family="symmetric", 
  			  control=loess.control(surface="direct"));

      if (visualize)
        lines(fit, lwd=2, col="red");
      
      # ii) Extrapolate the 
      G0 <- 0;
      R0 <- predict(fit, newdata=G0);

      offsetR <- c(G0=G0, R0=R0);
      if (visualize)
        points(x=offset[1], y=offset[2], col="red", lwd=2);

      #  i) Fit a loess curve through the data.
      fit <- loess(G ~ R, family="symmetric", 
  			  control=loess.control(surface="direct"));

      if (visualize)
        lines(x=fit$fitted, y=fit$x, lwd=2, col="green");
      
      # ii) Extrapolate the 
      R0 <- 0;
      G0 <- predict(fit, newdata=R0);

      offsetG <- c(G0=G0, R0=R0);
      if (visualize)
        points(x=offset[2], y=offset[1], col="green", lwd=2);

      offset <- rbind(offsetR, offsetG);
    }
    else if (method == "loess-lm") {
      #  i) Fit a loess curve through the data.
      fit <- loess(R ~ G, family="symmetric");

      if (visualize)
        lines(fit, lwd=2, col="red");
      
      g <- seq(0, 2e3);
      r <- predict(fit, newdata=g);
      fit <- lm(r ~ g);

      offsetR <- c(G0=0, R0=coef(fit)[1]);
      if (visualize)
        points(x=offset[1], y=offset[2], col="red", lwd=2);

      #  i) Fit a loess curve through the data.
      fit <- loess(G ~ R, family="symmetric");

      if (visualize) {
        x <- fit$fitted;
        y <- fit$x;
        o <- order(x);
        lines(x=x[o], y=y[o], lwd=2, col="green");
      }
      
      r <- seq(0, 2e3);
      g <- predict(fit, newdata=r);
      fit <- lm(g ~ r);

      offsetG <- c(G0=coef(fit)[1], R0=0);
      if (visualize)
        points(x=offset[2], y=offset[1], col="green", lwd=2);

      offset <- rbind(offsetR, offsetG);
    }

    res[[slide]] <- offset;
  } # for (slide...)

  res;
}, private=TRUE)  # findChannelBiasDifferences()


############################################################################
# HISTORY:
# 2006-02-08
# o Rdoc bug fix: Did not load SMA 'mouse.data' in as.MAData() of RGData.
# o Rd bug fix: Replaced section 'equals' with 'examples'.
# 2005-07-19
# o Replaced all path="" arguments to path=NULL.
# 2005-06-11
# o Small optimization of local log.na() in getM() and getA().
# 2005-03-02
# o Harmonized method arguments with generic functions.
# 2005-02-09
# o Added getChannelNames().
# 2005-02-02
# o Removed the recommendation to think of microarray data as (A,M).
# o Added getLogRatios() and getLogIntensities() and their Rdoc comments.
#   getM() and getA() are aliases for these.
# o Removed deprecated dyeSwap().
# 2004-02-17
# o Added getWithinChannelPairs() dated 2003-11-09.
# 2003-10-13
# o The PCA methods for findChannelBiasDifferences() does now only run
#   R vs G and not G vs R since they give the same result.
# 2003-09-18
# o Added private findChannelBiasDifferences().
# 2003-07-28
# o Added shift().
# 2003-06-15
# o BUG FIX: plotSpatial() did only work for "known" fields. User added 
#   fields did not plot.
# 2003-04-12
# o BUG FIX: getM() and getA() gave an error because SpotSlideData() did
#   not exist. Should be SpotSlideArray().
# o The SpotSlideArray objects returned by getM() and getA() do now also
#   include a reference to the Layout object.
# 2003-01-08
# o For convenience, the default values for arguments in plotXY() of RGData
#   has been changed; ylog=xlog and ylim=xlim.
# 2002-12-10
# o BUG FIX: as.RawData() did not work for data with only one slide.
# 2002-12-01
# o dyeSwap() for RawData was misstakenly redefined for the MAData class.
# 2002-10-24
# o Added normalizeQuantile().
# 2002-10-23
# o Added virtual fields M and A via the getM() and getA() method approach.
# 2002-10-14
# o Renamed dyeSwap() to swapDyes().
# 2002-05-29
# * RGData constructor now accept vectors for R and G by calling as.matrix(). 
# 2002-05-25
# * as.MAData() now sets M and A to NA if either R < 0 or G < 0.
#   Requested by Lei Jiang.
# 2002-05-05
# * BUG FIX: as.MAData() lost the slide names (colnames).
# 2002-05-04
# * BUG FIX: highlight() would not work since argument 'what' was not set.
# 2002-04-21
# * Added trial version of normalizeGenewise().
# 2002-04-20
# * Added a trial version of read(). Most of the job is done in support
#   functions in the MicroarrayData class.
# 2002-02-29
# * BUG FIX:  Sometimes as.MAData() was producing "integer overflow in R*G".
#   This was due to that R and G might be integers. Fixed by forcing the
#   multiplication to be done in double precision! Problem first reported
#   by lj22@u.washington.edu.
# 2002-02-26
# * Removed append(), as.data.frame(), get() etc since they are now
#   generic in the MicroarrayData class.
# * Rewrote code to make use of setMethodS3().
# * Updated the Rdoc's.
# 2002-01-24
# * Renamed all get() to extract().
# 2001-08-08
# * Updates dyeSwap with the argument 'slides' and wrote its Rdoc comments.
# 2001-08-07
# * Updated append to also add the layout if it is not already set.
# 2001-08-01
# * Added the field spot in as.data.frame().
# 2001-07-14
# * Update equals to include Layout comparison too. Added Rdoc for equals().
# 2001-07-11
# * Updated some of the Rdoc comments.
# 2001-07-02
# * Improved getColors().
# 2001-07-01
# * Removed plotSpatial(). Better to transform to MAData and plot there.
# 2001-06-28
# * Added a lot more Rdoc comments.
# * Rename argument 'indices' in highlight() to 'include' for consistancy
#   with other plot functions.
# 2001-05-14
# * Added getInternalReferences() for improving gco() performance.
# 2001-04-11
# * Added getColors().
# 2001-04-04
# * Added set/getXLabel() methods.
# * Made all cex=par("usr") by default.
# * Added highlight() to RGData also. First implemented in MAData.
# 2001-03-10
# * Created.
############################################################################
