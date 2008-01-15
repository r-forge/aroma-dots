#########################################################################/**
# @RdocClass MAData
#
# @title "The MAData class"
#
# \description{
#  @classhierarchy
#
#  Creates a new \code{MAData} object. 
#  The philosophy behind this data structure is to think about the data in the form of log ratios (M) and log intensites (A) for the spot signals. 
#  This is in contrast to the idea of the \link{RGData} structure, which thinks about the data as the signal in one channel (R) versus the signal in the other channel (G). 
# }
#
# @synopsis
#
# \arguments{
#   \item{M,A}{A NxM @matrix containing (base-2) log-ratios and
#    log-intensities, respectively.}
#   \item{layout}{A @see "Layout" object specifying the spot layout of the
#    slides in this data set.}
#   \item{extras}{Private argument. Do not use.}
# }
#
# \section{Fields and Methods}{
#  \bold{Fields}
#  \tabular{rll}{
#    \tab \code{M} \tab The log of the ratios between R and G, i.e. \eqn{\log_2{\frac{R}{G}}}{log2(R/G)} where R and G is the signal for the R and the G channel. \cr
#    \tab \code{A} \tab The log of the intensities of R and G, i.e. \eqn{\frac{1}{2}\log_2{R{\cdot}G}}{1/2*log2(R*G)} where R and G is the signal for the R and the G channel. \cr 
#  }
#
#  @allmethods
# }
#
# \details{
#   The mapping between M and A, and R and G is a one-to-one function, i.e.
#   you can go back and forth without loosing any information.
#   Given the signal R and G for the R and the G channels you get the
#   M and the A values by:
#   \deqn{
#     M = \log_2\frac{R}{G},\quad 
#     A = \log_2\sqrt{R{\cdot}G} = \frac{1}{2}\log_2 R{\cdot}G,
#   }{
#     M = log2(R/G), A = log2(sqrt(R*G)) = 1/2*log2(R*G),
#   }
#   and going back to the R and the G by:
#   \deqn{
#     R = \sqrt{2^{2A+M}},\quad G = \sqrt{2^{2A-M}}
#   }{
#     R = sqrt(2^(2A+M)), G = sqrt(2^(2A-M))
#   }
#
#   Comments: M is a memnonic for Minus and A is for Add since 
#   \eqn{M = \log_2{R}-\log_2{G}} and \eqn{A = 1/2*(\log_2{R}+\log_2{G})}.
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
#   # The option 'dataset' is used to annotate plots.
#   options(dataset="sma:mouse.data")
#
#   # Create a raw data object from the preexisting example data in
#   # the sma package.
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#
#   # Get the signal (here by default non-background corrected)
#   ma <- getSignal(raw, bgSubtract=TRUE)
#
#   # Transform (M,A) into (R,G)
#   rg <- as.RGData(ma)
#
#   # Transform back from (R,G) to (M,A)
#   ma2 <- as.MAData(rg);
#
#   # Check that the tranformation a one-to-one function
#   print(equals(ma, ma2))   # TRUE
#
#   layout(matrix(1:4, ncol=2, byrow=TRUE))
#   # Plot the R vs G with a fitted (lowess) line for slide 2.
#   plot(rg, slide=2); lowessCurve(rg)
#
#   # And the similar for M vs A.
#   plot(ma, slide=2); lowessCurve(ma)
#
#   # Plot a spatial representation of the M's.
#   plotSpatial(ma, slide=2)
#
#   # Make a boxplot of the print-tip groups.
#   boxplot(ma, groupBy="printtip", slide=2)
# }
#*/#########################################################################
setConstructorS3("MAData", function(M=NULL, A=NULL, layout=NULL, extras=list()) {
  if (!is.null(M)) M <- SpotSlideArray(M, layout=layout);
  if (!is.null(A)) A <- SpotSlideArray(A, layout=layout);
  
  this <- extend(MicroarrayData(layout=layout, extras=extras), "MAData",
    M      = M,
    A      = A,
    .cache = list()
  );
  this$.fieldNames <- c("M", "A");

  # Sets the default labels:
  setLabel(this, "M", expression(M==log[2](R/G)));
  setLabel(this, "A", expression(A==1/2*log[2](R*G)));

  this;
})



#########################################################################/**
# @RdocMethod swapDyes
#
# @title "Dye swap one or many slides"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \value{
#   Returns itself.
# }
#
# \arguments{
#   \item{slides}{A @vector of slides to be dye swapped. If @NULL, all
#   slides are considered.}
# }
#
# \examples{
#   # The option 'dataset' is used to annotate plots.
#   options(dataset="sma:MouseArray")
#
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#
#   ma <- getSignal(raw, bgSubtract=TRUE)
#   rg <- as.RGData(ma)
#
#   # Dye swap every other slide.
#   swapDyes(ma, slides=c(4,5,6))
#
#   layout(matrix(1:6, nrow=2, ncol=3, byrow=TRUE));
#   for (k in 1:6)
#     plot(ma, slide=k)
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("swapDyes", "MAData", function(this, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  
  this$M[,slides] <- -this$M[,slides];  

  clearCache(this);
  
  invisible(this);
})  # swapDyes()



setMethodS3("as.character", "MAData", function(this) {
  s <- paste(data.class(this), ": ", sep="");
  s <- paste(sep="",s,"M ",   com.braju.sma.dimStr(this$M));
  s <- paste(sep="",s,", A ", com.braju.sma.dimStr(this$A));
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


############################################################################
############################################################################
## 
##  DATA STRUCTURAL METHODS
## 
############################################################################
############################################################################

setMethodS3("as.MAData", "ANY", function(object, ...) {
  warning("Argument 'object' is of unknown type, but will try anyway.");
  MAData(obj, ...)
})


setMethodS3("as.MAData", "MAData", function(this, slides=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);

  M <- this$M[,slides];
  A <- this$A[,slides];

  layout <- getLayout(this);
  extras <- this$.extras;
  ma <- MAData(M=M, A=A, layout=layout, extras=extras)
  setSlideNames(ma, names=getSlideNames(this, slides=slides));
  ma;
})  # as.MAData()




#########################################################################/**
# @RdocMethod as.RGData
#
# @title "Transform MA format into RG format"
#
# \description{
#  Transform from log ratios and intensities to the red and green intensities.
#  Background fields are set to zero.
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{Subset of slides to be returned. If @NULL, all slides
#   are returned.}
# }
#
# \value{
#   Returns a @see "RGData" object.
# }
#
# \examples{
#   # The option 'dataset' is used to annotate plots.
#   options(dataset="sma:MouseArray")
#
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   rg <- getSignal(raw, bgSubtract=TRUE)
#   ma <- as.MAData(rg)
#   rg2 <- as.RGData(ma)
#   equal <- equals(rg, rg2)
# }
#
# \seealso{
#   @see "RGData.as.MAData".
#   @seeclass
# }
#
# @author
#*/#########################################################################
setMethodS3("as.RGData", "MAData", function(this, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  M <- this$M[,slides];
  A <- this$A[,slides];

  R <- as.matrix(sqrt(2^(2*A+M)));
  G <- as.matrix(sqrt(2^(2*A-M)));

  rg <- RGData(R=R, G=G, layout=getLayout(this), this$.extras)
  setSlideNames(rg, names=getSlideNames(this, slides=slides));
  rg;
})  # as.RGData()


setMethodS3("getR", "MAData", function(this) {
  M <- this$M;
  A <- this$A;
  SpotSlideArray(sqrt(2^(2*A+M)), layout=getLayout(this));
}, protected=TRUE)  # getR()

setMethodS3("getG", "MAData", function(this) {
  M <- this$M;
  A <- this$A;
  SpotSlideArray(sqrt(2^(2*A-M)), layout=getLayout(this));
}, protected=TRUE)  # getG()


setMethodS3("as.RawData", "MAData", function(this, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  rg <- as.RGData(this);
  R <- as.matrix(rg$R[,slides]);
  G <- as.matrix(rg$G[,slides]);

  Rb <- matrix(0, nrow=nrow(R), ncol=ncol(R));
  Gb <- Rb;

  raw <- RawData(R=R, G=G, Rb=Rb, Gb=Gb, layout=getLayout(this),
                                                extras=this$.extras)
  setSlideNames(raw, names=getSlideNames(this, slides=slides));
  raw;
})  # as.RawData()



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
# @title "Generates colors for each of the specified spots"
#
# \description{
#  Generates colors for each of the specified spots, which can be passed as
#  a \code{col} argument in most plot functions.
# }
#
# @synopsis
#
# \arguments{
#  \item{slide}{Specifies for which slide the colors should be generated for.}
#  \item{include}{The spot indices to be plotted. If @NULL all spots are considered.}
#  \item{exclude}{The spot indices \emph{not} to be plotted. If @NULL no spots are excluded (from all or from the one specified by \code{include}).}
#  \item{palette}{The palette to be used. Currently, only the \code{"redgreen"} palette is available.}
#  \item{M.range}{The range of "normal" M values. All M values outside this range will be cutoff.}
#  \item{A.range}{The range of "normal" A values. All A values outside this range will be cutoff.}
# }
#
# \value{Returns a @vector of colors.}
#
# @author
#
# \examples{
#   # The option 'dataset' is used to annotate plots.
#   options(dataset="sma:MouseArray")
#
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw, bgSubtract=TRUE)
#   col <- getColors(ma, slide=4)
# }
#
# \seealso{
#   @see "R.graphics::Colors".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getColors", "MAData", function(this, what=c("A","M"), slide=1, include=NULL, exclude=NULL, palette=NULL, range="auto", ...) {
  if (is.null(palette)) palette <- "auto";
  if (identical(what, c("M","A")))
    what <- c("A", "M");

  if (identical(palette, "auto")) {
    if (identical(what, c("A","M")))
      palette <- "redgreen"
    else if (all(is.element(c("R","G"), what)))
      palette <- "redgreen"
    else if (identical(what, "M"))
      palette <- "redgreen"
    else
      palette <- "grayscale";
  }

  if (identical(range, "auto")) {
    if (identical(what, c("A","M")))
      range <- matrix(c(0,16, -2,2), nrow=2)
    else if (identical(what, "A"))
      range <- c(0,16)
    else if (identical(what, "M"))
      range <- c(-2,2)
    else if (all(is.element(c("R","G"), what)))
      range <- matrix(c(0,16, 0,16), nrow=2)
    else
      range <- NULL;
  }

  include <- this$getInclude(include=include, exclude=exclude, slide=slide);

  x <- NULL;
  if (palette == "redgreen") {
    if (is.null(dim(range)) || dim(range) != c(2,2))
      range <- matrix(c(0,16, -2,2), nrow=2);
    if (all(is.character(what))) {
      for (field in c("A", "M"))
        x <- cbind(x, this[[field]][include,slide]);  # Don't forget slide here!
      if (identical(what, "M")) 
        x[,1] <- 13;
    } else {
      x <- what;
    }
    x[is.na(x)] <- 0;
    dim.range <- matrix(c(0,1, Colors$GREEN.HUE,Colors$RED.HUE), nrow=2);
    colors <- Colors$getHSV(x, x.range=range, 
                                   dim=c("v", "h"), dim.range=dim.range);
  } else if (palette == "grayscale") {
    if (all(is.character(what))) {
      what <- what[1];
      x <- this[[what]][include,slide];
    } else {
      x <- what;
    }
    if (is.null(range))
      range <- range(x, na.rm=TRUE)
    else
      range <- as.matrix(range)[,1];
    x[is.na(x)] <- 0;
    colors <- Colors$getGray(x, x.range=range);
  } else {
    throw("Unknown palette: ", palette);
  }
  colors;
})  # getColors()



setMethodS3("getDenseSpots", "MAData", function(this, x="A", y="M", nclass=10) {
  x <- this[[x]];
  y <- this[[y]];

  # Move origin to "lower left corner" of data set.
  x <- x-min(x, na.rm=TRUE);
  y <- y-min(y, na.rm=TRUE);

  # Normalize to [0,1] in (x,y).
  x <- x/max(x, na.rm=TRUE);
  y <- y/max(y, na.rm=TRUE);

  x2y2 <- x^2+y^2;
  n <- length(x2y2);

  h <- hist(x2y2, nclass=nclass, plot=FALSE);
  dense.bins <- (h$counts > n/nclass);
  dense.spots <- (1:n)[abs(x2y2-h$mids[dense.bins]) < 1/nclass];
  dense.spots;
})  # getDenseSpots()





#########################################################################/**
# @RdocMethod plotMvsM
#
# @title "Plots the log-ratios for one slide versus another"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{slides}{@vector of two slide indices to be plotted against 
#   each other.}
#  \item{pch}{The plot symbol.}
#  \item{xlim, ylim}{The visible range on the x and the y axis.}
#  \item{xlab, ylab}{The labels of the x and the y axis.}
#  \item{...}{Other arguments accepted by @see "graphics::plot".}
# }
#
# @author
#
# @examples "../incl/MAData.plotMvsM.Rex"
#
# \seealso{
#   @seemethod "pointsMvsM".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("plotMvsM", "MAData", function(this, slides=NULL, include=NULL, exclude=NULL, pch="auto", xlim=NULL, ylim=xlim, xlab=NULL, ylab=NULL, ..., style=NULL) {
  returnValue <- NULL;

  cmd <- NULL;
  if (!is.null(style) && is.element(style, c("points", "highlight", "text", "lowessline"))) {
    cmd <- style;
    lastPlot <- Device$getPlotParameters();
    if (is.null(slides))
      slides <- lastPlot$slides;
  }

  slides <- validateArgumentSlides(this, slides=slides);
  if (length(slides) != 2) {
    throw("Argument 'slides' should specify exactly two slides: ", 
                                        paste(slides, collapse=", "));
  }

  # Needed? /HB 2003-12-30
  setView(this, MicroarrayArray$DEFAULT.VIEW);

  include <- which(getInclude(this, include, exclude, slide=slides[1]));

  if (is.null(pch)) pch <- par("pch");

  slideNames <- getSlideNames(this, slides=slides);
  if (is.null(slideNames))
    slideNames <- as.character(slides);
  if (is.null(xlab))
    xlab <- substitute(M^{(s)}, list=list(s=slideNames[1]));
  if (is.null(ylab))
    ylab <- substitute(M^{(s)}, list=list(s=slideNames[2]));

  Mx <- this$M[include,slides[1]]; 
  My <- this$M[include,slides[2]]; 

  Device$setPlotParameters(object=this, fcn="plotMvsM", slides=slides);

  if (is.null(cmd)) {  
    if (length(pch) == 1 && pch == "auto")
      pch <- 176;
    plot(Mx,My, pch=pch, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...);
    Device$putTimestamp();
    Device$putDataset();
    Device$putAuthor();
  } else if (cmd == "points") {
    if (length(pch) == 1 && pch == "auto")
      pch <- 176;
    points(Mx,My, pch=pch, ...);
  } else if (cmd == "highlight") {
    if (length(pch) == 1 && pch == "auto")
      pch <- 20;
    points(Mx,My, pch=pch, ...);
  } else {
    throw("Not supported: ", cmd);
  }
})


setMethodS3("plot", "MAData", function(x, what="MvsA", xlim=NULL, ylim=NULL, ...) {
  # To please R CMD check...
  this <- x;

  if (identical(what, "MvsA") && is.null(ylim)) {
    # Assert equal scale on the M and the A axis. This will keep R and G ortogonal
    # in the M vs A plot.

    # By making 'ylim' into an expression it can be evaluated later.
    ylim <- expression(c(-1,1)*(xlim[2]-xlim[1]));
  }
  if (is.null(xlim) && is.null(ylim)) {
    plot.MicroarrayData(this, what=what, ...);
  } else if (is.null(xlim)) {
    plot.MicroarrayData(this, what=what, ylim=ylim, ...);
  } else if (is.null(ylim)) {
    plot.MicroarrayData(this, what=what, xlim=xlim, ...);
  } else {
    plot.MicroarrayData(this, what=what, xlim=xlim, ylim=ylim, ...);
  }
})

setMethodS3("plotXY", "MAData", function(this, what=c("A","M"), ...) {
  plotXY.MicroarrayData(this, what=what, ...);
})

setMethodS3("plotSpatial", "MAData", function(this, what=c("A","M"), ...) {
  plotSpatial.MicroarrayData(this, what=what, ...);
})

setMethodS3("plotDensity", "MAData", function(this, what="M", ...) {
  NextMethod("plotDensity", this, what=what, ...);
})

setMethodS3("boxplot", "MAData", function(x, what="M", ...) {
  # To please R CMD check...
  this <- x;

  boxplot.MicroarrayData(this, what=what, ...);
})

setMethodS3("hist", "MAData", function(x, what="M", ...) {
  # To please R CMD check...
  this <- x;

  hist.MicroarrayData(this, what=what, ...);
})

setMethodS3("qqnorm", "MAData", function(y, slide=1, include=NULL, exclude=NULL, ...) {
  # To please R CMD check...
  this <- y;

  include <- this$getInclude(include=include, exclude=exclude, slide=slide);
  M <- this$M[include];

  qqnorm(M, ...);
  qqline(M);
})


setMethodS3("plotPrintorder", "MAData", function(this, what="M", ...) {
  plotPrintorder.MicroarrayData(this, what=what, ...);
})


setMethodS3("plotDiporder", "MAData", function(this, what="M", ...) {
  plotDiporder.MicroarrayData(this, what=what, ...);
})


hideme01 <- function() {
      if (!is.null(extra.type)) {
      # Extra on top of that
  	if(extra.type == "t")
  	    txt(A,M, labels, pch=pch.ex, cex=cex.ex, col=col.ex, ...);
  	if(extra.type == "p")
  	    points(A,M, pch=pch.ex, cex=cex.ex, col=col.ex, ...);
    
  	if(extra.type == "tci") { # Text labels
  	    M.ex <- M[topN];
  	    A.ex <- A[topN];
  	    if (length(cex.ex) > 1) 
  	      cex.ex <- cex.ex[topN];
  	    if (length(col.ex) > 1)
  	      col.ex <- col.ex[topN];
  	    if (length(pch.ex) > 1) 
  	      pch.ex <- pch.ex[topN];
    
  	    # Place text labels on top and below actual value.
  	    pos.ex <- sign(M.ex)+2;
  	    
  	    text(A.ex, M.ex, labels=labels[topN], cex=cex.ex, col=col.ex, pos=pos.ex, ...);
  	    points(A.ex, M.ex, pch.ex=pch.ex, cex=cex.ex, ...);
  	    str(cex.ex);
  	}
    
  	if(extra.type == "pci")
  	    plot.confband.points(A,M, crit1,crit2, nclass,
  				      pch=pch.ex, cex=cex.ex, col=col.ex, ...);
  	if(extra.type == "lci") 
  	    plot.confband.lines(A,M, crit1,crit2, nclass,
  				      pch=pch.ex, cex=cex.ex, col=col.ex, ...);
      } # if (!is.null(extra.type))
    
      # Plot the lowess line...
      ind <- !(is.na(A) | is.na(M) | is.infinite(A) | is.infinite(M))
      lowess.line <- lowess(A[ind], M[ind], f=0.3);
      lines(lowess.line, lty=1, col="black")
      lines(lowess.line, lty=3, col="white")
}





############################################################################
############################################################################
## 
##  STATISTICAL METHODS
## 
############################################################################
############################################################################

setMethodS3("range", "MAData", function(this, what="M", slide=NULL, na.rm=TRUE, inf.rm=TRUE, ...) {
  if (is.element(what, getFields(this))) {
    X <- this[[what]];
    if (!is.null(slide)) {
      if (!is.matrix(X))
        throw("The field that 'what' refers to is not a applicable field: ", what);
      X <- X[,slide];
    }
    if (inf.rm) X[is.infinite(X)] <- NA;
    range(X, na.rm=na.rm);
  } else
    range.MicroarrayData(this, what=what, na.rm=na.rm, inf.rm=inf.rm, ...);
})  # range()



setMethodS3("getRange", "MAData", function(this, what=c("M", "A"), slides=NULL) {
  warning("getRange() is deprecated. Use range() instead!!!");
  slides <- validateArgumentSlides(this, slides=slides);

  if (!is.element(what, getFields(this)))
    throw("The argument 'what' does not specify a data field: ", what);

  X <- this[[what]];
  if (!is.matrix(X))
    throw("The field that 'what' refers to is not a applicable field: ", what);

  X <- as.matrix(X[,slides]);
  X.range <- apply(X, MARGIN=2, FUN=range);
  c(min(X.range[1,]),max(X.range[2,]));
}, protected=TRUE)   # getRange()




#########################################################################/**
# @RdocMethod mean
#
# @title "Average Mean for microarray data"
#
# \description{
#   Computes the average mean across all slides for each spot.
# }
#
# @synopsis
#
# \arguments{
#  \item{inf.rm}{a logical value indicating whether @Inf values should be
#         stripped before the computation proceeds.}
#  \item{na.rm}{a logical value indicating whether @NA values should be
#         stripped before the computation proceeds.}
#  \item{...}{other arguments to \code{mean}, e.g. \code{trim}.}
# }
#
# \value{Returns a @matrix with the columns \code{meanM} and \code{meanA}.}
#
# \examples{
#   # The option 'dataset' is used to annotate plots.
#   options(dataset="sma:MouseArray")
#
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw, bgSubtract=TRUE)
#   mean <- mean(ma)
# }
#
# @author
#
# \seealso{
#   @seemethod "var".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("mean", "MAData", function(x, inf.rm=TRUE, na.rm=TRUE, ...) {
  # To please R CMD check...
  this <- x;

  fields <- c("M", "A");
  res <- matrix(0, nrow=nbrOfSpots(this), ncol=length(fields));
  colnames(res) <- paste("v", fields, sep="");
  for (k in 1:length(fields)) {
    x <- this[[fields[k]]];
    if (inf.rm) x[is.infinite(x)] <- NA;
    res[,k] <- apply(x, MARGIN=1, FUN=mean, na.rm=na.rm, ...);
  }
  res;
})  # mean()




#########################################################################/**
# @RdocMethod var
#
# @title "Variance for microarray data"
#
# \description{
#   Computes the variance across all slides for each spot.
# }
#
# @synopsis
#
# \arguments{
#  \item{inf.rm}{a logical value indicating whether @Inf values should be
#         stripped before the computation proceeds.}
#  \item{na.rm}{a logical value indicating whether @NA values should be
#         stripped before the computation proceeds.}
#  \item{...}{other arguments to \code{var}, e.g. \code{trim}.}
# }
#
# \value{Returns a @matrix with the columns \code{varM} and \code{varA}.}
#
# \examples{
#   # The option 'dataset' is used to annotate plots.
#   options(dataset="sma:MouseArray")
#
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw, bgSubtract=TRUE)
#   var <- ma$var()
# }
#
# @author
#
# \seealso{
#   @seemethod "mean".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("var", "MAData", function(this, inf.rm=TRUE, na.rm=TRUE, ...) {
  fields <- c("M", "A");
  res <- matrix(0, nrow=nbrOfSpots(this), ncol=length(fields));
  colnames(res) <- paste("v", fields, sep="");
  for (k in 1:length(fields)) {
    x <- this[[fields[k]]];
    if (inf.rm) x[is.infinite(x)] <- NA;
    res[,k] <- apply(x, MARGIN=1, FUN=var, na.rm=na.rm, ...);
  }
  res;
})  # var()









setMethodS3("getHistogram", "MAData", function(this, what="M", slide=1, ...) {
  if (what == "M") {
    data <- this$M[,slide];
  } else if (what == "A") {
    data <- this$A[,slide];
  } else
    throw("Argument 'what' must be either \"M\" or \"A\": ", what);

  hist(unclass(data), plot=FALSE, ...);
})  # getHistogram()



#########################################################################/**
# @RdocMethod topSpots
#
# @title "Gets the top spots"
#
# @synopsis
#
# \description{
#  Gets (approximately) the \code{n} spots with the highest M values, if 
#  \code{n} greater or equal to one, and the 100*\code{n} percentage if 
#  \code{n} is less than one.
# }
#
# @author
#
# \examples{
#   # The option 'dataset' is used to annotate plots.
#   options(dataset="sma:MouseArray")
#
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw, bgSubtract=TRUE)
#   normalizeWithinSlide(ma, "s")
#   top100  <- topSpots(ma, 100);   # The top 100 spots
#   top0.01 <- topSpots(ma, 0.01);  # The top 1 percent of the spots
# }
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("topSpots", "MAData", function(this, n=50, nbins=1, slide=1, 
            include=NULL, exclude=NULL, limits=c("both", "upper", "lower"), 
                                                             X="A", Y="M") {
  if (is.null(X) && nbins > 1)
    throw("Argument 'X' must be given if 'nbins' is greater than one.");

  ratio <- (n < 1);

  if (limits == "both")
    mult <-2
  else
    mult <-1;

  if (!ratio && n < mult*nbins)
    throw("Number of spots (n) wanted must be at least as many as the number of bins (nbins), since at least two spots are choosen from each bin.");

  include <- this$getInclude(include=include, exclude=exclude, slide=slide);

  y <- this[[Y]][include,slide];
  if (is.null(X))
    x <- 1:length(y)
  else
    x <- this[[X]][include,slide];

  # Removes values that are infinite
  x[is.infinite(x)] <- NA;
  y[is.infinite(y)] <- NA;

  # Divide spots into nbins equally sized bins over x
  xbin <- quantile(x, probs=seq(0, nbins, 1)/nbins, na.rm=TRUE);

  if (!ratio) {
    # Nbr of spots in each bin to pick
    nbin <- n/nbins;
    # Half on top and half on bottom... (depends on 'limits')
    nbin2 <- nbin/mult;
  }

  spots <- rep(FALSE, length(include));

  # For each bin, find the topn top spots
  for(i in 1:nbins) {
    # Choose among y's that are in the current bin...
    ybin <- y;
    bin <- (xbin[i] <= x) & (x < xbin[i+1]);
    ybin[!bin] <- NA;
    if (ratio) {
      topn <- n;
    } else {
      # Percentage of spots wanted for this bin.
      len <- length(ybin[!is.na(ybin)]);
      topn <- nbin2 / len;
    } 

    if (limits == "both") {
      prob0 <- topn;
      prob1 <- 1-topn;
    } else if (limits == "upper") {
      prob0 <- 0;
      prob1 <- 1-topn;
    } else {
      prob0 <- topn;
      prob1 <- 1;
    }

    # Select those spots in current bin that fullfills the criteria...
    cutoff <- quantile(ybin, probs=c(prob0, prob1), na.rm=TRUE);
    vals <- ((ybin < cutoff[1]) | (ybin > cutoff[2]));
    spots[vals] <- TRUE;
  }
  
  include[spots];
})  # topSpots()



setMethodS3("plotSpatial3d", "MAData", function(this, field="M", col=NULL, log=NULL, ...) {
  if (is.null(col)) {
    if (identical(field, "M"))
      what <- c("A", "M")
    else
      what <- field;
    col <- getColors(this, what=what, log=log);
  }
  plotSpatial3d.MicroarrayData(this, field=field, col=col, log=log, ...);
})  # plotSpatial3d()


setMethodS3("plot3d", "MAData", function(this, ...) {
  plotSpatial3d(this, ...);
})






#########################################################################/**
# @RdocMethod getGeneVariability
#
# @title "Gets the genewise variability"
#
# @synopsis
#
# \arguments{
#   \item{robust}{If @TRUE the median absolute deviation (MAD) of
#     the residuals will be calculated, otherwise the sample standard
#     deviation will be calculated.}
#   \item{force}{If @FALSE and if cached gene variability values 
#     exists they will be used, otherwise the gene variability will be
#     (re-)calculated.}
#   \item{slides}{The slides which should be included in the calculations.
#     If @NULL, all slides are included.}
# }
#
# \description{
#   Calculates the genewise variability. What spots belongs to which genes
#   is defined by the layout of the slides, where all slides are assumed to
#   have the same layout. See class @see "Layout" for more
#   information about genes. For speed improvement, the gene variability 
#   will be cached for future queries. To override the cache, use 
#   \code{force=TRUE}.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw, bgSubtract=TRUE)
#   ma.norm <- clone(ma)
#   normalizeWithinSlide(ma.norm, method="s")
#   normalizeAcrossSlides(ma.norm)
#   var <- getGeneVariability(ma)
#   var.norm <- getGeneVariability(ma.norm)
#   var.diff <- var - var.norm
#
#   # Statistics
#   print(summary(var))
#   print(summary(var.norm))
#   print(summary(var.diff))
#   cat("Number of 'improved' genes :", sum(var.diff > 0, na.rm=TRUE), "\n")
#   cat("Number of 'worsened' genes :", sum(var.diff < 0, na.rm=TRUE), "\n")
#   cat("Number of 'unchanged' genes:", sum(var.diff == 0, na.rm=TRUE), "\n")
#
#   # Plots
#   xlim <- quantile(c(var, var.norm), probs=c(0,0.999))
#   subplots(4)
#   hist(var, nclass=100, xlim=xlim, cex.main=0.7,
#             main="Gene variability before normalization");
#   hist(var.norm, nclass=100, xlim=xlim, cex.main=0.7,
#                  main="Gene variability after normalization");
#   hist(var.diff, nclass=50, cex.main=0.7, main="Improvements");
# }
#
# @author
#
# \seealso{
#   To calculate the mean (or any other quantile) of the genewise 
#   variabilities see @seemethod "getMOR".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getGeneVariability", "MAData", function(this, robust=TRUE, force=FALSE, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  geneVariability <- getCache(this, "geneVariability", force=force);
  if (!is.null(geneVariability))
    return(geneVariability);
  
  # Compare this with the SE in T = mean(M) / SE(M)
  layout <- getLayout(this);
  geneGroups <- getGeneGroups(layout);
  genes <- getSpots(geneGroups);
  ngenes <- nbrOfGroups(geneGroups);
  # as.matrix() will retrieve the matrix of the SpotSlideArray,
  # which will speed up the calculations below a lot, because
  # there is quite a bit of overhead in "[.SpotSlideArray".
  data  <- this$M[,slides];

  variability <- rep(NA, ngenes);
  if (robust == TRUE) {
    # Robust standardization
    for (l in 1:ngenes) {
      spots <- genes[[l]];
      x <- data[spots,];
      ok <- !is.na(x);
      x <- x - median(x[ok]);
      variability[l] <- 1.4826*median(abs(x[ok]));
    }
  } else {
    # Non-robust standardization
    for (l in 1:ngenes) {
      spots <- genes[[l]];
      x <- data[spots,];
      ok <- !is.na(x);
      variability[l] <- sqrt(var(x[ok]));  # Standard Deviation
    }
  } # if (robust == ...)

  names(variability) <- names(genes);
  attr(variability, "df") <- as.numeric(getSizes(geneGroups));
  class(variability) <- "GeneVariability";

  # Set the cache
  setCache(this, "geneVariability", variability);
})  # getGeneVariability()


setMethodS3("getSpotVariability", "MAData", function(this, robust=TRUE, force=FALSE, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  variability <- getCache(this, "spotVariability", force=force);
  if (!is.null(variability))
    return(variability);
  
  # Compare this with the SE in T = mean(M) / SE(M)
  nspots <- nbrOfSpots(this);
  # as.matrix() will retrieve the matrix of the SpotSlideArray,
  # which will speed up the calculations below a lot, because
  # there is quite a bit of overhead in "[.SpotSlideArray".
  data  <- this$M[,slides];

  variability <- rep(NA, nspots);
  if (robust == TRUE) {
    # Robust standardization
    variability <- 1.4826 * apply(data, MARGIN=1, FUN=function(x) {
      ok <- !is.na(x);
      x <- x - median(x[ok]);
      median(abs(x[ok]));
    })
  } else {
    # Non-robust standardization
    variability <- apply(data, MARGIN=1, FUN=function(x) {
      ok <- !is.na(x);
      sqrt(var(x[ok]));
    })
  } # if (robust == ...)

  names(variability) <- 1:nspots;
  attr(variability, "df") <- length(slides);
  class(variability) <- "SpotVariability";

  # Set the cache
  setCache(this, "spotVariability", as.matrix(variability));
})  # getSpotVariability()



#########################################################################/**
# @RdocMethod getMOR
#
# @title "Gets the Measure of Reproducibility"
#
# @synopsis
#
# \arguments{
#   \item{robust}{If @TRUE the median absolute deviation (MAD) will
#     be used for calculating the genewise residuals, otherwise the sample
#     standard deviation will be used.}
#   \item{probs}{The quantiles of the genewise variabilities returned. If
#     \code{-0.5} (or @NULL) the mean is returned.}
#   \item{force}{If @FALSE and if cached gene variability values exists
#     they will be used, otherwise the gene variability will be
#     (re-)calculated.}
#   \item{slides}{The slides which should be included in the calculations.
#     If @NULL, all slides are included.}
# }
#
# \description{
#   Calculates the Measure of Reproducibility (MOR) [1], which is by default
#   the (scalar) mean value of all genewise robust variability measures
#   given by \code{getGeneVariability()}, either the standard deviation or
#   the median absolute deviation (MAD). Any quantile(s) of these can be
#   returned by setting argument \code{probs}.
# }
#
# \value{
#   Returns a @vector of the quantiles of the genewise variabilities asked
#   for. By default the mean value of all genewise variabilities (MOR) is
#   returned.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw, bgSubtract=TRUE)
#   ma.norm <- clone(ma)
#   normalizeWithinSlide(ma.norm, method="s")
#   normalizeAcrossSlides(ma.norm)
#   mor <- getMOR(ma)
#   mor.norm <- getMOR(ma.norm)
#
#   cat("MOR for non-normalized data:", mor, "\n")
#   cat("MOR for normalized data:", mor.norm, "\n")
# }
#
# @author
#
# \references{
#  [1] Henrik Bengtsson, Plate Effects in cDNA microarray data, 
#      Matemathical Statistics, Centre for Matematical Sciences,
#      Lund University, Sweden. Manuscript, 2002.
# }
#
# \seealso{
#   @seemethod "getVariability"
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getMOR", "MAData", function(this, robust=TRUE, probs=NULL, force=FALSE, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (is.null(probs)) probs <- -0.5;
  meanIdx <- (probs == -0.5);
  ps <- probs[!meanIdx];
  if (length(ps) > 0 && (range(ps)[1] < 0 || range(ps)[2] > 1))
    throw("Argument 'probs' is out of range: ", probs, collapse=", ");

  d <- getGeneVariability(this, robust=robust, force=force, slides=slides);

  mor <- rep(NA, length(probs));
  mor[meanIdx] <- mean(d, na.rm=TRUE);
  names(mor)[meanIdx] <- "mean";
  
  if (length(ps) > 0) {
    q <- quantile(d, probs=ps, na.rm=TRUE);
    mor[!meanIdx] <- q;
    names(mor)[!meanIdx] <- names(q);
  }

  mor;
})  # getMOR()




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
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw, bgSubtract=TRUE)
#
#   subplots(4)
#   xlim <- c(4,16); ylim <- c(-3,3);
#   plot(ma, xlim=xlim, ylim=ylim)
#   min1 <- function(x) { min(x)-1 }   # Shift to signal one (not zero!)
#   shift(ma, R=min, G=min)
#   plot(ma, xlim=xlim, ylim=ylim)
#   shift(ma, M=median)
#   plot(ma, xlim=xlim, ylim=ylim)
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("shift", "MAData", function(this, M=NULL, A=NULL, R=NULL, G=NULL, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  foo <- function(x, X) {
    if (is.function(X)) {
      x <- x - X(x[is.finite(x)]);
    } else {
      x <- x + X;
    }
  }

  for (slide in slides) {
    # i) shift M and A
    if (!is.null(M))
      this$M[,slide] <- foo(this$M[,slide], M);
    if (!is.null(A))
      this$A[,slide] <- foo(this$A[,slide], A);

    # ii) shift R and G
    if (!is.null(R) || !is.null(G)) {
      m <- this$M[,slide];
      a <- this$A[,slide];
      r <- sqrt(2^(2 * a + m));
      g <- sqrt(2^(2 * a - m));
      if (!is.null(R))
  	r <- foo(r, R);
      if (!is.null(G))
  	g <- foo(g, G);
      m <- log(r/g, base=2);
      a <- 1/2*log(r*g, base=2);
      rm(r,g);
      this$M[,slide] <- m;
      this$A[,slide] <- a;
    }
  } # for (slide in slides)

  clearCache(this);
}) # shift()

setMethodS3("shiftEqualRG", "MAData", function(this, a=NULL, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (is.character(a)) {
    if (a == "max") {
      # Shift so that the the smallest signal is equal to one.
      minR <- apply(getR(this), MARGIN=2, FUN=min, na.rm=TRUE);
      minG <- apply(getG(this), MARGIN=2, FUN=min, na.rm=TRUE);
      min <- minG;
      min[minR < minG] <- minR[minR < minG];
      aMax <- -min+1;
      for (slide in slides)
        shift(this, R=aMax, G=aMax, slides=slide);
    } else {
      throw("Unknown shift method: ", a);
    }
  } else {
    shift(this, R=a, G=a, slides=slides);
  }

  clearCache(this);
}) # shiftEqualRG()


#########################################################################/**
# @RdocMethod loessQuantile
#
# @title "Estimates the (25\%,50\%,75\%) quantile M(A) curves using loess"
#
# @synopsis
#
# \arguments{
#   \item{slides}{The slides which should be included in the calculations.
#     If @NULL, all slides are included.}
#   \item{subset}{Indices of spots to be used for the estimates. Selecting
#     a random subset may speed up the estimation and still give the same
#     result. If @NULL, all spots are included.}
#   \item{method}{If \code{"each"}, estimates of the quantiles are given
#     for each of the slides. If \code{"average"}, an estimate of the 
#     quantiles is given for all spots as if they came from the same slide.}
#   \item{factor}{@numeric factor to multiply to the upper and lower 
#     quantiles. If the data is Gaussian distributed the default factor
#     \code{1.4826} makes the upper and lower quantile to correspond to
#     one standard deviation.}
#   \item{family}{Family to be used by @see "modreg::loess". 
#     By default, a robust family is used.}
#   \item{...}{Other arguments passed to \code{loess()}.}
# }
#
# \description{
#   @get "title". 
#   This method is used by @seemethod "getAdjustedSpotVariability".
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#
#   # Create four sets of slides where slide 2 and 4 are shifted R=G=a
#   a <- 2^11
#   ma <- list()
#   ma[[1]] <- getSignal(raw, bgSubtract=TRUE)
#   ma[[2]] <- clone(ma[[1]])
#   shift(ma[[2]], R=a, G=a)
#   ma[[3]] <- clone(ma[[1]])
#   normalizeWithinSlide(ma[[3]], method="p")
#   ma[[4]] <- clone(ma[[3]])
#   shift(ma[[4]], R=a, G=a)
#
#   setMethodS3("lines", "LoessQuantile", function(object, ...) {
#     lines(object$x, object$centre, ...);
#     lines(object$x, object$above, ...);
#     lines(object$x, object$below, ...);
#   })
#
#   # Four plots
#   subplots(4)
#   Alim <- Mlim <- NA
#   for (k in 1:length(ma)) {
#     Alim <- range(c(Alim, ma[[k]]$A), na.rm=TRUE)
#     Mlim <- range(c(Mlim, ma[[k]]$M), na.rm=TRUE)
#   }
#   
#   for (k in 1:length(ma)) {
#     plot(ma[[k]], slide=1, xlim=Alim, ylim=Mlim, col=k+1)
#     fit <- loessQuantile(ma[[k]], slide=1)[[1]]
#     lines(fit, lwd=2)
#   }
# }
#
# @author
#
# \seealso{
#   @seemethod "getAdjustedSpotVariability".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("loessQuantile", "MAData", function(this, slides=NULL, subset=NULL, method=c("each", "average"), factor=1.4826, family="symmetric", ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  method <- match.arg(method);

  loessQuantileOne <- function(A,M, ...) {
    predictA <- seq(from=min(A), to=max(A), length.out=100);
  
    centre <- loess(M ~ A, family=family, ...);
  
    # Identify spots that are above the centre line.
    Mcentre <- predict(centre, newdata=A);
    up <- (M-Mcentre > 0);
  
    # Keep only a few data points to save memory
    centre <- predict(centre, newdata=predictA);
  
    M <- M - Mcentre;
    above <- loess(M ~ A, subset=up, family=family, ...);
    above <- predict(above, newdata=predictA);
  #  above <- centre + factor*(above-centre);
    above <- centre + factor*above;
  
    below <- loess(M ~ A, subset=!up, family=family, ...);
    below <- predict(below, newdata=predictA);
  #  below <- centre + factor*(below-centre);
    below <- centre + factor*below;
  
    scale <- (above - below);
  
    res <- list(x=predictA, centre=centre, above=above, below=below, scale=scale);
    attr(res, "class") <- "LoessQuantile";
    res;
  } # loessQuantileOne()

  if (is.null(subset)) {
    subset <- rep(TRUE, length.out=nbrOfSpots(this));
  } else {
    if (is.numeric(subset)) {
      subset0 <- rep(FALSE, length.out=nbrOfSpots(this));
      subset0[subset] <- TRUE;
      subset <- subset0;
    } else if (!is.logical(subset)) {
      throw("Argument 'subset' must be logical or numeric: ", mode(subset));
    }
  }

  if (method == "each") {
    res <- list();
    for (slide in slides) {
      A <- this$A[,slide];
      M <- this$M[,slide];
      ok <- (is.finite(A) & is.finite(M));
      ok <- ok & subset;
      A <- A[ok];
      M <- M[ok];
      gc();
      res[[slide]] <- loessQuantileOne(A=A, M=M, ...);
    } # for (slide ...)
    res;
  } else {
    A <- as.vector(this$A[,slides]);
    M <- as.vector(this$M[,slides]);
    ok <- (is.finite(A) & is.finite(M));
    ok <- ok & subset;
    A <- A[ok];
    M <- M[ok];
    gc();
    loessQuantileOne(A=A, M=M);
  }
}, protected=TRUE) # loessQuantile()




#########################################################################/**
# @RdocMethod getAdjustedSpotVariability
#
# @title "Gets the spotwise intensity-adjusted variability of replicate slides"
#
# @synopsis
#
# \arguments{
#   \item{robust}{If @TRUE the median absolute deviation (MAD) of
#     the residuals will be calculated, otherwise the sample standard
#     deviation will be calculated.}
#   \item{force}{If @FALSE and if cached gene variability values 
#     exists they will be used, otherwise the gene variability will be
#     (re-)calculated.}
#   \item{slides}{The slides which should be included in the calculations.
#     If @NULL, all slides are included.}
# }
#
# \description{
#   @get "title". Within-slide replicates are considered to be independent
#   of each other.
# }
#
# @examples "../incl/MAData.getAdjustedSpotVariability.Rex"
#
# @author
#
# \seealso{
#   See also @seemethod "getSpotVariability" for non-intensity dependent
#   scale adjustment.
#   variabilities see @seemethod "getMOR".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getAdjustedSpotVariability", "MAData", function(this, robust=TRUE, force=FALSE, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  loessQuantile <- getCache(this, "loessQuantile", force=force);
  if (is.null(loessQuantile)) {
    subset <- sample(1:nbrOfSpots(this), 1000);
    loessQuantile <- loessQuantile(this, slides=slides, subset=subset, method="average");
    setCache(this, "loessQuantile", loessQuantile);
  }

  d <- getSpotVariability(this, robust=robust, force=force, slides=slides);
  medA <- apply(this$A[,slides], MARGIN=1, FUN=median, na.rm=TRUE);
  scale <- loessQuantile$scale;
  x <- loessQuantile$x;
  fit <- loess(scale ~ x);
  ok <- is.finite(medA);
  s <- rep(NA, length.out=length(medA));
  s[ok] <- predict(fit, newdata=medA[ok]);
  w <- 1/s;
#  w <- w / mean(w, na.rm=TRUE);
  list(dw=d*w, d=d, w=w);
}, protected=TRUE) # getAdjustedSpotVariability()



#########################################################################/**
# @RdocMethod getMOR2003a
#
# @title "Gets the Intensity Adjusted Measure of Reproducibility"
#
# @synopsis
#
# \arguments{
#   \item{robust}{If @TRUE the median absolute deviation (MAD) will
#     be used for calculating the genewise residuals, otherwise the sample
#     standard deviation will be used.}
#   \item{probs}{The quantiles of the genewise variabilities returned. If
#     \code{-0.5} (or @NULL) the mean is returned.}
#   \item{force}{If @FALSE and if cached gene variability values exists
#     they will be used, otherwise the gene variability will be
#     (re-)calculated.}
#   \item{slides}{The slides which should be included in the calculations.
#     If @NULL, all slides are included.}
# }
#
# \description{
#   Calculates the Adjusted Measure of Reproducibility (MOR2003a), which is
#   by default the (scalar) mean value of all genewise robust variability 
#   measures given by \code{getGeneVariability()}, either the standard
#   deviation or the median absolute deviation (MAD),  weighted by the 
#   inverse of the estimated standard deviation of the log-ratios at the 
#   given log-intensities. Any quantile(s) of these can be returned by
#   setting argument \code{probs}.
# }
#
# \value{
#   Returns a @vector of the quantiles of the genewise variabilities asked
#   for. By default the mean value of all genewise variabilities (MOR2003a) 
#   is returned.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   # Keep only slides of treatment 1
#   mouse.data <- lapply(mouse.data, FUN=function(x) x[,4:6])
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw, bgSubtract=TRUE)
#   ma.norm <- clone(ma)
#   normalizeWithinSlide(ma.norm, method="s")
#   normalizeAcrossSlides(ma.norm)
#   mor <- getMOR2003a(ma)
#   mor.norm <- getMOR2003a(ma.norm)
#
#   cat("MOR2003a for non-normalized data:", mor, "\n")
#   cat("MOR2003a for normalized data:", mor.norm, "\n")
# }
#
# @author
#
# \seealso{
#   @seemethod "getAdjustedSpotVariability"
#   @seemethod "getMOR"
#   @seemethod "getSpotVariability"
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getMOR2003a", "MAData", function(this, robust=TRUE, probs=NULL, force=FALSE, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(probs)) probs <- -0.5;
  meanIdx <- (probs == -0.5);
  ps <- probs[!meanIdx];
  if (length(ps) > 0 && (range(ps)[1] < 0 || range(ps)[2] > 1))
    throw("Argument 'probs' is out of range: ", probs, collapse=", ");

  d <- getAdjustedSpotVariability(this, robust=robust, force=force, slides=slides);
  d <- d$dw;

  mor <- rep(NA, length(probs));
  mor[meanIdx] <- mean(d, na.rm=TRUE);
  names(mor)[meanIdx] <- "mean";
  if (length(ps) > 0) {
    q <- quantile(d, probs=ps, na.rm=TRUE);
    mor[!meanIdx] <- q;
    names(mor)[!meanIdx] <- names(q);
  }

  mor;
}, protected=TRUE) # getMOR2003a()





############################################################################
# HISTORY:
# 2006-06-13
# o BUG FIX: From R v2.3.1 graphics::image.default() does not accept 
#   xlim=NULL nor ylim=NULL.  Had to update plot.MAData() accordingly.
# 2006-02-08
# o Rd bug fix: Replaced section 'amauments' with 'arguments'.
# 2005-02-02
# o Removed deprecated dyeSwap().
# 2004-02-17
# o Added plotDensity() dated 2003-11-10.
# 2004-01-08
# o For plot() in MAData the ylim will now follow the xlim such that 
#   R and G, which are rotate 45 degrees in the M vs A plot, are still
#   ortogonal to each other. This is recommended because the logarithmic 
#   scale fools you enough anyway.
# 2003-12-30
# o Added plotMvsM() date 2003-11-04.
# 2003-10-13
# o BUG FIX: Assigned cache in getGeneVariability() with wrong value (NULL).
#   Lead to setCache(obj, name, NULL), which in turn gave an error, which
#   it shouldn't! Fixed too.
# 2003-09-24
# o BUG FIX: getAdjustedSpotVariability() could not handle cases where 
#   slides contain log-intensities with value NA.
# 2003-09-17
# o Added protected methods getAdjustedSpotVariability(), getMOR2003a() and 
#   the support method loessQuantile(). Wrote example code for them too.
# 2003-07-28
# o Added shift().
# 2003-07-07
# o as.BMAData(): Added links and recommendations about the limma package
#   instead.
# 2003-04-12
# o getR() and getG() now returns a SpotSlideArray with a reference to 
#   the Layout object too.
# 2003-04-08
# o BUG FIX: as.TMAData() gave "Error in var.default(Mi) : `x' is empty" if
#   the matrix M contains a row with all NA's. Thanks Valtteri Wirta at
#   the Royal Institute of Technology, Stockholm for this bug report.
# 2002-12-10
# o BUG FIX: as.RawData() did not work for data with only one slide.
# 2002-12-06
# o Now getGeneVariability() returns a n-by-1 matrix.
# 2002-11-13
# o Added argument 'slides' to getGeneVariability() and getMOR().
# o Make sure to call unclass(this$M) whenever you just need the matrix. 
#   Any "[" operator is much faster on the matrix.
# 2002-11-12
# o All fields are now of class SpotSlideArray.
# 2002-10-23
# o Added virtual fields R and G via the getR() and getG() method approach.
# 2002-10-22
# o Package updated to make use of new R.oo.
# 2002-10-14
# o Renamed dyeSwap() to swapDyes().
# 2002-10-09
# o Update as.TMAData() to be make use of getGeneGroups(), which will always
#   guarantee that t-statistics are calculated genewise, not spotwise.
# 2002-09-30
# o Made getGeneResiduals() deprecated and private.
# o Renames getGeneVariance() to getGeneVariability(), because it was the
#   MAD or the standard deviation that was calculated and not the variance.
# o Added Rdoc comments to getMOR() and made it public.
# 2002-06-30
# * Splitted MAData.R into MAData.R and MAData.NORM.R.
# * BUG FIX: getMOR() gave an error if probs==-0.5 ("mean") only.
# 2002-06-24
# * BUG FIX: mean() was broken.
# * Made getGeneDiscrepancies() obsolete.
# * Added argument 'weights' to normalizeWithinSlide() on request from Jon
#   McAuliffe. Currently only 0-1 weights or FALSE-TRUE weights are
#   supported, since lowess() does not support weights. If loess() would
#   be used weights could be [0,1].
# 2002-05-30
# * Added trial versions of getGeneVariance() and getMOR().
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
