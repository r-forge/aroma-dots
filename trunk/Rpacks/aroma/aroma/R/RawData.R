#########################################################################/**
# @RdocClass RawData
#
# @title "The RawData class"
#
# \description{
#  @classhierarchy
#
#  Creates a new \code{RawData} object.
# }
#
# @synopsis
#
# \arguments{
#   \item{R,G}{A NxM @matrix containing (non-logged) foreground signals of
#    the red (green) channel, where N is the number of spots on each slide
#    and M is the number of slides in this data set.}
#   \item{Rb,Gb}{A NxM @matrix containing (non-logged) background signals of
#    the red (green) channel, where N is the number of spots on each slide
#    and M is the number of slides in this data set.}
#   \item{layout}{A @see "Layout" object specifying the spot layout of the
#    slides in this data set.}
#   \item{extras}{Private argument. Do not use.}
# }
#
# \section{Fields and Methods}{
#  \tabular{rll}{
#    \tab \code{R} \tab The foreground for channel R (non-logged). \cr
#    \tab \code{G} \tab The foreground for channel G (non-logged). \cr
#    \tab \code{Rb} \tab The background for channel Rb (non-logged). \cr
#    \tab \code{Gb} \tab The background for channel Gb (non-logged). \cr
#  }
#
#  @allmethods
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
# \examples{\dontrun{For an example see help(MAData).}}
#*/#########################################################################
setConstructorS3("RawData", function(R=NULL, G=NULL, Rb=NULL, Gb=NULL, 
                                                 layout=NULL, extras=list()) {
  # This is to support old sma style too.
  if (is.list(R)) {
    if ( !any(is.na(match(c("R", "G", "Rb", "Gb"), names(R)))) ) {
      l <- R;
      R  <- l$R;
      G  <- l$G;
      Rb <- l$Rb;
      Gb <- l$Gb;
    }
  }

  if (!is.null(R))   R <- SpotSlideArray(R);
  if (!is.null(G))   G <- SpotSlideArray(G);
  if (!is.null(Rb)) Rb <- SpotSlideArray(Rb);
  if (!is.null(Gb)) Gb <- SpotSlideArray(Gb);

  this <- extend(MicroarrayData(layout=layout, extras=extras), "RawData",
    R  = R,
    G  = G,
    Rb = Rb,
    Gb = Gb
  );
  this$.fieldNames <- c("R", "G", "Rb", "Gb");
  
  # Sets the default labels:
  setLabel(this, "R", expression(R));
  setLabel(this, "G", expression(G));
  setLabel(this, "Rb", expression(Rb));
  setLabel(this, "Gb", expression(Gb));

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
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#
#   # Dye swap every other slide.
#   swapDyes(raw, slides=c(4,5,6))
#
#   layout(matrix(1:6, nrow=2, ncol=3, byrow=TRUE));
#   for (k in 1:6)
#     plot(raw, "RvsG", slide=k)
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("swapDyes", "RawData", function(this, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  tmp <- this$R[,slides];
  this$R[,slides] <- this$G[,slides];  
  this$G[,slides] <- tmp;
  tmp <- this$Rb[,slides];
  this$Rb[,slides] <- this$Gb[,slides];  
  this$Gb[,slides] <- tmp;

  clearCache(this); 
  
  invisible(this);
})



setMethodS3("as.character", "RawData", function(this) {
  s <- "RawData: ";
  s <- paste(sep="",s,"R ",    com.braju.sma.dimStr(this$R));
  s <- paste(sep="",s,", G ",  com.braju.sma.dimStr(this$G));
  s <- paste(sep="",s,", Rb ", com.braju.sma.dimStr(this$Rb));
  s <- paste(sep="",s,", Gb ", com.braju.sma.dimStr(this$Gb));
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

setMethodS3("as.RawData", "RawData", function(this, slides=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);

  R  <- this$R[,slides];
  G  <- this$G[,slides];
  Rb <- this$Rb[,slides];
  Gb <- this$Gb[,slides];

  layout <- getLayout(this);
  extras <- this$.extras;
  RawData(R=R, G=G, Rb=Rb, Gb=Gb, layout=layout, extras=extras)
})


setMethodS3("as.RawData", "ANY", function(this, ...) {
  RawData(this, ...)
})




#########################################################################/**
# @RdocMethod getSignal
#
# @title "Gets the red and the green signal"
#
# \description{
#   @get "title" as a MAData object.
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{Subset of slides to be returned. If @NULL, all slides
#   are returned.}
#   \item{bgSubstract}{If @TRUE, the background is subtracted from the
#         foreground, before the transformation is performed.
#         This argument was previously named \code{bg.subtract}, which still
#         works, but is deprecated and will be removed in a future version.}
# }
#
# \value{
#   Returns a @see "MAData" object.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw, bgSubtract=FALSE);
# }
#
# @author
#
# \seealso{
#   @seemethod "getBackground",
#   @seemethod "getForeground".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getSignal", "RawData", function(this, bgSubtract=TRUE, slides=NULL, bg.subtract=TRUE) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (missing(bgSubtract) && !missing(bg.subtract)) {
    bgSubtract <- bg.subtract;
    warning("Argument 'bg.subtract' in getSignal() (class RawData) is deprecated. Please use 'bgSubtract' instead.");
  }

  R <- this$R[,slides];
  G <- this$G[,slides];
  if (bgSubtract) {
    R <- R - this$Rb[,slides];
    G <- G - this$Gb[,slides];
  }

  layout <- getLayout(this);
  extras <- this$.extras;
  rg <- RGData(R=R, G=G, layout=layout, extras=extras)
  ma <- as.MAData(rg);
})



#########################################################################/**
# @RdocMethod getForeground
#
# @title "Gets the foreground signal"
#
# \description{
#  @get "title" as a MAData object.
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
#   Returns the foreground signal as a @see "MAData" object.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   fg <- getForeground(raw)
#   plot(fg)
# }
#
# @author
# 
# \seealso{
#   @seemethod "getBackground",
#   @seemethod "getSignal".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getForeground", "RawData", function(this, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  R <- this$R[,slides];
  G <- this$G[,slides];

  layout <- getLayout(this);
  extras <- this$.extras;
  rg <- RGData(R=R, G=G, layout=layout, extras=extras)
  as.MAData(rg);
})



#########################################################################/**
# @RdocMethod getBackground
#
# @title "Gets the background signal"
#
# \description{
#  @get "title" as a MAData object.
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
#   Returns the background signal as a @see "MAData" object.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   bg <- getBackground(raw)
#   plot(bg)
# }
#
# @author
# 
# \seealso{
#   @seemethod "getForeground",
#   @seemethod "getSignal".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getBackground", "RawData", function(this, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  R <- this$Rb[,slides];
  G <- this$Gb[,slides];

  layout <- getLayout(this);
  extras <- this$.extras;
  rg <- RGData(R=R, G=G, layout=layout, extras=extras)
  as.MAData(rg);
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
#   Returns the background signal as a @see "RawData" object.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   print(raw)
#   rawGPairs <- getWithinChannelPairs(raw)
#   print(rawGPairs)
# }
#
# @author
# 
# \seealso{
#   See also \code{getSlidePairs()} in the @see "MicroarrayData@ class,
#   which is used internally.
#   @seeclass
# }
############################################################################
setMethodS3("getWithinChannelPairs", "RawData", function(this, channel, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  pairs <- getSlidePairs(this, slides=slides);

  channelBg <- paste(channel, "b", sep="");

  X1 <- this[[channel]][,pairs[1,]];
  X2 <- this[[channel]][,pairs[2,]];
  Xb1 <- this[[channelBg]][,pairs[1,]];
  Xb2 <- this[[channelBg]][,pairs[2,]];

  RawData(R=X1,G=X2,Rb=Xb1,Gb=Xb2, layout=getLayout(this));
}) # getWithinChannelPairs()





############################################################################
############################################################################
## 
##  PLOTTING & GRAPHICAL METHODS
## 
############################################################################
############################################################################

setMethodS3("plotXY", "RawData", function(this, what=c("Rb", "R"), xlog=2, ylog=xlog, xlim=NULL, ylim=xlim, ...) {
  plotXY.MicroarrayData(this, what, xlog=xlog, ylog=ylog, xlim=xlim, ylim=ylim, ...);
})


setMethodS3("plotSpatial", "RawData", function(this, what="R", col="auto", ...) {
  if (!is.null(col)) {
    if (col == "auto") {
      whats <- c("R", "Rb", "G", "Gb");
      cols <- c("redscale", "redscale", "greenscale", "greenscale");
      col <- cols[what == whats];
    } else {
      col <- "grayscale";
    }
  }
  plotSpatial.MicroarrayData(this, what=what, col=col, ...);
})


setMethodS3("boxplot", "RawData", function(x, what="R", ...) {
  # To please R CMD check...
  this <- x;

  boxplot.MicroarrayData(this, what=what, ...);
})


setMethodS3("plot", "RawData", function(x, what="RvsRb", ...) {
  # To please R CMD check...
  this <- x;

  plot.MicroarrayData(this, what=what, ...);
})



#########################################################################/**
# @RdocMethod getColors
# 
# @title "Generates colors for each of the specified spots"
#
# \description{
#  @get "title", which can be passed to
#  the \code{col} argument in most plot functions. However, note that most
#  of the plot functions in this package automatically recognize palettes
#  if they are passed to the \code{col} argument.
# }
#
# @synopsis
#
# \arguments{
#  \item{what}{The fields to be used.}
#  \item{slide}{Specifies for which slide the colors should be generated for.}
#  \item{include}{The indices of the spots that should be included. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are considered.
#   If @NULL all spots are considered.}
#  \item{exclude}{The indices of the spots that should be excluded. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are excluded.}
#   \item{bgSubstract}{If @TRUE, the background is subtracted from the
#         foreground, before the colors are generated.}
#   \item{palette}{The palette to be used. If @NULL or \code{"auto"} 
#         the palette that is most applicable to the value of \code{what}
#         will be used.}
#   \item{range}{Cut-off range for data values. The format is a matrix with
#     two rows and \code{length(what)} number of columns.}
# }
#
# \value{Returns a @vector of colors.}
#
# @author
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   col <- getColors(raw, what="G", slide=4);
# }
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getColors", "RawData", function(this, what=c("R","G"), slide=1, include=NULL, exclude=NULL, palette=NULL, bgSubtract=FALSE, bg.subtract=bgSubtract, range=matrix(c(0,16,0,16), nrow=2), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  log.na <- function (x, ...) {
    ifelse(x > 0, log(x, ...), NA)
  }

  if (missing(bgSubtract) && !missing(bg.subtract)) {
    bgSubtract <- bg.subtract;
    warning("Argument 'bg.subtract' in getColors() (class RawData) is deprecated. Please use 'bgSubtract' instead.");
  }

  if (is.null(palette)) palette <- "auto";

  bgs <- na.omit(match(c("Rb", "Gb"), what));
  if (length(what) == 2 && length(bgs) != 1) {
    if (length(bgs) == 0)
      rg <- getForeground(this)
    else
      rg <- getBackground(this)
    ma <- as.MAData(rg);
    getColors(ma, what=c("A","M"), slide=slide, include=include, exclude=exclude, palette=palette);
  } else {
    include <- getInclude(this, include=include, exclude=exclude, slide=slide);

    if (length(bgs) == 0) bgs <- 1;
    bgfield <- what[bgs];
    bg <- this[[bgfield]][include, slide];

    if (length(what) == 1) {
      sg <- bg;
      channel <- substring(bgfield,1,1);
    } else {
      sgfield <- what[-bgs][1];
      if (!is.element(sgfield, c("R", "G"))) {
  	whatStr <- paste(what, collapse=", ");
  	throw("Unknown value of argument what (\"", whatStr, "\").");
      }
      sg <- this[[sgfield]][include, slide];
      if (bgSubtract) sg <- sg-bg;
      channel <- sgfield;
    }

    x <- log.na(sg,2);
    x[is.na(x)] <- 0;

    palettes <- c("grayscale", "auto", "redgreen", "redscale", "greenscale");
    idx <- match(palette, palettes);
    if (is.null(idx))
      throw("Unknown palette (\"", palette, "\").");
    if (idx == 1) {
      Colors$getGray(x, x.range=range);
    } else {
      dimensions <- c("", channel, channel, "R", "G");
      dim <- dimensions[idx];
      Colors$getRGB(x, x.range=range, dim=dim);
    }
  }
})



############################################################################
############################################################################
## 
##  STATISTICAL METHODS
## 
############################################################################
############################################################################

setMethodS3("range", "RawData", function(this, what="R", slide=NULL, na.rm=TRUE, inf.rm=FALSE, ...) {
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


setMethodS3("read", "RawData", function(this, filename, path=NULL, layout=NULL, verbose=FALSE) {
  fields <- c("R", "G", "Rb", "Gb");
  res <- MicroarrayData$readToList(filename, path=path,
                                   reqFields=fields, verbose=verbose);
  
  # Create the MAData object.
  RawData(R=res$R, G=res$G, Rb=res$Rb, Gb=res$Gb, layout=layout)
}, static=TRUE, trial=TRUE);



  # This is the signal-to-noise ratios as defined by Wang et al. [1]
  #
  #  [1] Quantitative quality control in microarray image processing
  #      and data acquisition, X. Wang, S. Ghosh and S-W. Guo,
  #      Nucleic Acids Research, 2001, Vol. 29, No. 15 e75, 2001,
  #      Oxford University Press.
  #      \url{http://nar.oupjournals.org/cgi/content/abstract/29/15/e75}.

setMethodS3("getSNR1", "RawData", function(this, log=NULL) {
  R <- this$R;
  G <- this$G;
  Rb <- this$Rb;
  Gb <- this$Gb;
  if (!is.null(log)) {
    R  <- log(R, base=log);
    G  <- log(G, base=log);
    Rb <- log(Rb, base=log);
    Gb <- log(Gb, base=log);
  }
  R.SNR <- R / (R + Rb);
  G.SNR <- G / (G + Gb);
  SNR <- sqrt(R.SNR * G.SNR);           # Geometric mean
  attr(SNR, "R") <- R.SNR;
  attr(SNR, "G") <- G.SNR;
  SNR;
})

setMethodS3("getSNR2", "RawData", function(this, log=NULL) {
  R <- this$R;
  G <- this$G;
  Rb <- this$Rb;
  Gb <- this$Gb;
  if (!is.null(log)) {
    R  <- log(R, base=log);
    G  <- log(G, base=log);
    Rb <- log(Rb, base=log);
    Gb <- log(Gb, base=log);
  }
  R.SNR <- R / Rb;
  G.SNR <- G / Gb;
  SNR <- sqrt(R.SNR * G.SNR);           # Geometric mean
  attr(SNR, "R") <- R.SNR;
  attr(SNR, "G") <- G.SNR;
  SNR;
})






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
#   \item{M,A,R,G,Rb,Gb}{A @numeric or @function specifying the shift to be 
#    applied to the log-ratios, the log-intensities, the red, and/or the 
#    green foreground or background signals.
#    If more than one of these are shifted at the same time, they are
#    shifted in the order \code{M}, \code{A}, \code{R}, \code{G}, 
#    \code{Rb} and \code{Gb}.
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
setMethodS3("shift", "RawData", function(this, M=NULL, A=NULL, R=NULL, G=NULL, Rb=NULL, Gb=NULL, slides=NULL) {
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

    # ii) shift R and G
    if (!is.null(R))
      this$R[,slide] <- foo(this$R[,slide], R);
    if (!is.null(G))
      this$G[,slide] <- foo(this$G[,slide], G);

    # iii) shift Rb and Gb
    if (!is.null(Rb))
      this$R[,slide] <- foo(this$Rb[,slide], Rb);
    if (!is.null(Gb))
      this$G[,slide] <- foo(this$Gb[,slide], Gb);
    }
  } # for (slide in slides)

  clearCache(this); 
})



############################################################################
# HISTORY:
# 2005-07-19
# o Replaced all path="" arguments to path=NULL.
# 2005-06-11
# o Made log.na() a local function of getColors().
# 2005-02-02
# o Removed deprecated dyeSwap().
# 2004-02-17
# o Added getWithinChannelPairs() dated 2003-11-09.
# 2003-07-28
# o Added shift().
# 2003-06-15
# o BUG FIX: plotSpatial() did only work for "known" fields. User added 
#   fields did not plot.
# 2003-03-20
# o Made bg.subtract in getSignal() deprecated and replaced by bgSubtract.
# 2003-01-08
# o For convenience, the default values for arguments in plotXY() of RawData
#   has been changed; ylog=xlog and ylim=xlim.
# 2002-12-01
# o dyeSwap() for RawData was misstakenly redefined for the MAData class.
# 2002-10-24
# o Added normalizeQuantile().
# 2002-10-14
# o Renamed dyeSwap() to swapDyes().
# o BUG FIX: Updated the Rdoc field section.
# 2002-05-04
# o BUG FIX: highlight() would not work since argument 'what' was not set.
# 2002-04-21
# * Added trial version of normalizeGenewise()
# 2002-04-20
# * Added getSNR().
# * Added a trial version of read(). Most of the job is done in support
#   functions in the MicroarrayData class.
# 2002-02-26
# * Updated the Rdoc's.
# * Remove append, as.data.frame, extract, nbrOfSpots, nbrOfSlides etc
#   since they are now in a generic form in the class MicroarrayData.
# * Updated the code to make use of setMethodS3().
# 2002-01-24
# * Renamed all get() to extract().
# 2001-11-12
# * Added fromResultsData.RawData() which takes a list of ResultsData 
#   objects, which could have been generated by a readAll.ResultsData 
#   method. This should simplify the reading of several files at once.
# * Added "RGData:" in as.character().
# 2001-08-09
# * Added normalizeWithinSlide.
# 2001-08-08
# * Updates dyeSwap with the argument 'slides' and wrote its Rdoc comments.
# * Rdoc bug: Update usage of getSignal.
# 2001-08-07
# * Updated append to also add the layout if it is not already set.
# 2001-08-01
# * Added the field spot in as.data.frame().
# 2001-07-14
# * Update equals to include Layout comparison too. Added Rdoc for equals().
# 2001-07-12
# * In getSignal() bg.subtract is now by default TRUE.
# 2001-07-08
# * Now getSignal(), getForeground(), getBackground() return a MAData object
#   instead of a RGData object. The latter was hardly never used, but one
#   always went on to conform it to a MAData object.
# 2001-07-02
# * Bug fix: Updated getColors(). More safe now and supports 
#   include/exclude, palette etc.
# 2001-07-01
# * Came up with a new highlight approach which passes the highlighting
#   problem back to the last function.
# * Updated plotSignalvsBackground() a lot.
# 2001-05-14
# * Added getInternalReferences() for improving gco() performance.
# 2001-04-05
# * Added getForeground().
# 2001-04-03
# * Implemented the plotSignalvsBackground() myself.
# 2001-03-27
# * Removed plot() that plotted R vs G. Use getSignal()$plot() or 
#   getBackground()$plot() instead.
# * Removed plot.sb(). Non-informative name!
# 2001-03-23
# * Created.
############################################################################

