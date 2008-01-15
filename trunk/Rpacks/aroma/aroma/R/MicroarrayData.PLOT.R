############################################################################
# MicroarrayData - Plot Functions
#
# Method:           Req #fields: View(s): Description:  
# ------------------------------------------------------------------------
# plot()            2            default  scatter plot
# [plotXY()]                     (ss,gs)
#
# plotSpatial()     1            ss       2D spatial plot
# plotSpatial3d()   1            ss       3D spatial plot
#
# plotDensity()     1            ss       empirical density distributions
#
# plotPrintorder()  1            ss       print-order plot
# plotDiporder()    1            ss       dip-order plot, i.e. average 
#                                         print-order plot for all p-tips.
#
# plotGene()        1            grs      plot within and between slides
#                                         replicated values for one gene
# plotReplicates()  1            grs      plot within and between slides
# [obsolete?]                             replicated values for one gene
#
# boxplot()         1            default  boxplot grouped by ???
#                                (ss,gs)
# hist()            1            default  histogram
#                                (ss,gs)
# qqnorm()          1            default  quantile vs normal quantile plot
#                                (ss,gs)  
#
# Special plot methods:
# highlight()       as last      as last  highlight data points in the last
#                                         generated plot.
# points()          as last      as last  adds points to the last generated
#                                         plot.
# text()            as last      as last  adds labels to data points in the 
#                                         last generated plot.
# lowessCurve()     as last      as last  adds a lowess curve to the last
#                                         generated plot.
#
# To do:
# legend()          as last      as last  adds a legend for the last plot.
# plotMvsA()
# plotYvsY()
# plotMvsM()
# plotMcvsAc()
#
# Legends: 
#   ss  - (spot,slide) view, 
#   gs  - (gene,slide) view, 
#   grs - (gene,replicate,slide) view
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#  Plot methods for LayoutGroups:
#  a) Several groups in one plot:
#   plotGroups()*
#
#  b) One group in one plot:
#   plotGroup()*
#
#  subplots()
#
#  createColors()
#  getColors()
#
# Auxillary (protected) methods:
#  plotXY()
#
# *) To be implemented.
############################################################################


setMethodS3("isFieldColorable", "MicroarrayData", function(this,  field) {
  is.element(field, c("R", "G", "Rb", "Gb", "Ravg", "Gavg", 
                      "M", "A", "Mavg", "Aavg", "T"));
})

setMethodS3("setView", "MicroarrayData", function(this, newView) {
  this$.view <- newView;
})

setMethodS3("getView", "MicroarrayData", function(this) {
  this$.view;
})




#########################################################################/**
# @set "class=MicroarrayData"
# @RdocMethod plot
#
# @title "Plots spatial representation of a microarray spot statistics"
#
# \description{
#  @get "title". Currently
#  the statistics can be either the log-ratios (default) or the intensities.
# }
#
# @synopsis
#
# \arguments{
#  \item{what}{What to plot. Common formats are \code{"YvsX"} and
#   \code{"spatialX"} where \code{X} and \code{Y} are variables in the
#   object.}
#  \item{...}{Common arguments accepted by most plot functions.
#   For more information see @seemethod "plotXY", @seemethod "plotSpatial",
#   @see "graphics::par", and @see "graphics::plot".}
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#   rg <- as.RGData(ma)
#
#   layout(matrix(1:6, nrow=2, byrow=TRUE))
#   plot(raw)            # Plots R vs Rb. Default is "RvsRb" and slide 1.
#   plot(raw, "GvsGb")   # Plots G vs Gb.
#   plot(rg)             # Plots R vs G. Default is "RvsG" and slide 1.
#   plot(rg, "GvsR")     # Plots G vs R.
#   plot(ma, slide=4)    # Plots M vs A for slide 4. Default is "MvsA".
#   plot(ma, "spatialM") # Plots spatial plot of M's.
# }
#
# \seealso{
#   @seemethod "plotXY" and @seemethod "plotSpatial", 
#   which are called by \code{plot}.
#   For highlightning and put labels on spots after plotting see 
#   @seemethod "highlight" and
#   @seemethod "text".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("plot", "MicroarrayData", function(x, what, ..., style=NULL) {
  # To please R CMD check...
  this <- x;

  if (!is.null(style) && is.element(style, c("points", "highlight", "text", "lowessline"))) {
    lastPlot <- Device$getPlotParameters();
    plotfcn <- get(lastPlot$fcn, mode="function");
    if (!is.null(what)) {
      searches <- c("vs", "spatial", "boxplot", "histogram", "printorder");
      match <- which(apply(as.matrix(searches), MARGIN=1, 
                           FUN=function(x) regexpr(x, what)) != -1);
      if (length(match) == 0)
        throw("Unknown value of argument 'what': ", what); 
      plotfcn <- plotfcns[[match]];
      what <- as.character(what);
      what <- rev(unlist(strsplit(what, searches[match])));
      what <- as.character(what);
      what <- unlist(strsplit(what, "&"));  # For supporting "histogramR&G".
      what <- what[nchar(what) != 0];
    }
  } else {
    searches <- c("vs", "spatial", "boxplot", "histogram", "printorder");
    plotfcns <- c(plotXY, plotSpatial, boxplot, hist, plotPrintorder);

    match <- which(apply(as.matrix(searches), MARGIN=1, 
                           FUN=function(x) regexpr(x, what)) != -1);
    if (length(match) == 0)
      throw("Unknown value of argument 'what': ", what);

    plotfcn <- plotfcns[[match]];
    what <- as.character(what);
    what <- rev(unlist(strsplit(what, searches[match])));
    what <- as.character(what);
    what <- unlist(strsplit(what, "&"));  # For supporting "histogramR&G".
    what <- what[nchar(what) != 0];
  }
  plotfcn(this, what=what, style=NULL, ...);
})





#########################################################################/**
# @RdocMethod plotSpatial3d
#
# @title "Plots a 3-dimensional spatial representation of a field"
#
# \description{
#  @get "title". 
# }
#
# @synopsis
#
# \arguments{
#  \item{field}{Name of the field to be plotted.}
#  \item{slide}{Slide number to be plotted.}
#  \item{pch}{Default value is 176 (small circle).}
#  \item{theta, phi}{Angles defining the viewing direction. 
#        \code{theta} gives the azimuthal direction and 
#        \code{phi} the colatitude.}
#  \item{xlab,ylab,zlab}{Labels for the x, the y, and the z axis.}
#  \item{log}{The base of logarithm to be used for the z dimension.
#             If @NULL, the logarithm is not calculated.}
#  \item{...}{Other arguments accepted by @see "R.basic::plot3d".}
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#   subplots(4)
#   plotSpatial3d(raw, "Rb", pch=".", col="red")
#   plotSpatial3d(raw, "Gb", pch=".", col="green")
#   plotSpatial3d(ma, "M", pch=".")
#   plotSpatial3d(ma, "A", pch=".")
# }
#
# \seealso{
#   @seemethod "plotXY" and @seemethod "plotSpatial", 
#   which are called by \code{plot}.
#   @see "R.basic::plot3d".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("plotSpatial3d", "MicroarrayData", function(this, field, slide=1, pch=176, phi=40, theta=-25, xlab="x", ylab="y", zlab=field, log=NULL, ...) {
  slide <- validateArgumentSlide(this, slide=slide);
  
  xy <- getSpotPosition(this, slide=slide);
  if (is.null(xy)) {
    layout <- getLayout(this);
    xy <- getPosition(layout);
    xy <- list(x=xy[,2], y=xy[,1]);
  }
  z <- this[[field]][,slide];
  if (!is.null(log))
    z <- log(z, base=log);
  plot3d(xy$x, -xy$y, z, xlab=xlab, ylab=ylab, zlab=zlab,
         pch=pch, phi=phi, theta=theta, ...);
}) # plotSpatial3d()



# Kept for backward compatibility
setMethodS3("plot3d", "MicroarrayData", function(this, ...) {
  warning("plot3d() in class MicroarrayData is deprecated.\n");
  plotSpatial3d(this, ...);
}, deprecated=TRUE);




#########################################################################/**
# @RdocMethod plotXY
#
# @title "Plots a scatter plot of two fields"
#
# \description{
#  @get "title".
#  It is recommended to use the \code{plot} function instead of calling
#  this method explicitly (see @seemethod "plot").
# }
#
# @synopsis
#
# \arguments{
#  \item{what}{What to plot. Any two field that can be retrieved by 
#   \code{extract} are accepted.}
#  \item{slide}{The slide to be plotted.}
#  \item{include}{The indices of the spots that should be included. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are considered.
#   If @NULL all spots are considered.}
#  \item{exclude}{The indices of the spots that should be excluded. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are excluded. If @NULL no spots are excluded.}
#  \item{gridwise}{If @TRUE a lowess line for each grid will be drawn,
#   otherwise only the global lowess line will be drawn.}
#  \item{xlog, ylog}{The logarithmic base to be used to take the
#   logarithm of the x and the y values, respectively. If @NULL, the
#   non-logarithmic values are plotted.
#   Note the difference from the definition in @see "graphics::par".}
#  \item{xlab, ylab}{The labels on the x and the y axis, respectively. If
#   @NULL, the fields default label will be used.}
#  \item{cex}{The scale factor to be used. If @NULL the system default
#   scale factor will be used.}
#  \item{col}{The color(s) to be used for the plotted spots, i.e. for the
#   spots \emph{after} inclusion and exclusion. If the value is
#   \code{"redgreen"} a red to green palette is used.}
#  \item{pch}{The dot style. Default value is \code{176}, which is a small
#   circle.}
#  \item{f}{The bandwidth for the lowess lines.}
#  \item{...}{Common arguments accepted by most plot functions.
#   For more information see @see "graphics::par" and @see "graphics::plot".}
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   plotXY(ma)                    # Plot M vs A (default)
#   plotXY(ma, what=c("M","A"))   # Plot A vs M.
# }
#
# @author
#
# \seealso{
#   @seemethod "plot".
#   @seemethod "plotSpatial".
#   @seemethod "highlight".
#   @seemethod "text".
#   @see "R.graphics::plotSymbols.Device".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("plotXY", "MicroarrayData", function(this, what, style=NULL, slide=NULL, include=NULL, exclude=NULL, xlog=NULL, ylog=NULL, xlab=NULL, ylab=NULL, cex=NULL, col="auto", pch="auto", gridwise=FALSE, f=0.3, lim=NULL, xlim=NULL, ylim=NULL, ...) {
  returnValue <- NULL;

  cmd <- NULL;
  if (!is.null(style) && is.element(style, c("points", "highlight", "text", "lowessline"))) {
    cmd <- style;
    lastPlot <- Device$getPlotParameters();
    if (is.null(what))
      what <- lastPlot$what;
    if (is.null(slide))
      slide <- lastPlot$slide;
    xlog <- lastPlot$xlog;
    ylog <- lastPlot$ylog;
  }

  if (length(what) != 2 || !is.character(what))
    throw("Argument 'what' must be a vector of two field names.");

  slide <- validateArgumentSlide(this, slide=slide);

  X <- what[1]; Y <- what[2];

  setView(this, MicroarrayArray$DEFAULT.VIEW);

  # The following row is a potential problem, because it *sorts* for instance
  # 'include' if it is a vector of spot indices. Then we have to sort 'col'
  # 'pch' etc too! I knew about this problem a long time ago, but forgot
  # about it again. Have to find a intuitive design for it! /HB 2004-02-23
  include <- which(getInclude(this, include, exclude, slide=slide));

  if (is.null(cex))   	   cex  <- par("cex");
  if (is.null(pch))   	   pch  <- par("pch");
  if (is.null(xlab))  	   xlab <- getLabel(this, X);
  if (is.null(ylab))  	   ylab <- getLabel(this, Y);
  if (is.expression(xlab)) xlab <- xlab[[1]];
  if (is.expression(ylab)) ylab <- ylab[[1]];

  # To make this generic method a little bit faster, first try to find the
  # data variables X and Y among the fields of the object, i.e. this[[X]]
  # and this[[Y]]. If they don't exists go an use the get() method instead.
  if (any(is.na(match(c(X,Y), names(this))))) {
    df <- extract(this, fields=c(X,Y), slide=slide);
    x.all <- df[,X];
    y.all <- df[,Y];
  } else {
    field <- this[[X]];

    field <- pushView(field, MicroarrayArray$DEFAULT.VIEW);
    x.all <- field[,slide];
    field <- popView(field);
    field <- this[[Y]];
    field <- pushView(field, MicroarrayArray$DEFAULT.VIEW);
    y.all <- field[,slide];
    field <- popView(field);
  }

  if (!is.null(xlog)) {
    x.all <- log(x.all, base=xlog);
    xlab <- substitute(log[YY](XX), list(XX=xlab, YY=xlog));
  }
  if (!is.null(ylog)) {
    y.all <- log(y.all, base=ylog);
    ylab <- substitute(log[YY](XX), list(XX=ylab, YY=ylog));
  }

  if (length(col) == 1) {
    if (!is.numeric(col) && substring(col,1,1) != "#" && !is.element(col, colors()))
      col <- getColors(this, what=what, slide=slide, include=include, palette=col);
  }

  x <- x.all[include];
  y <- y.all[include];

  if (is.null(xlim)) {
    if (is.null(lim)) {
      xlim <- range(x[is.finite(x)]);
    } else {
      xlim <- lim;
    }
  }

  if (is.null(ylim)) {
    if (is.null(lim)) {
      ylim <- range(y[is.finite(y)]);
    } else {
      ylim <- lim;
    }
  }

  if (is.expression(xlim))
     xlim <- eval(xlim);
  if (is.expression(ylim))
     ylim <- eval(ylim);

  Device$setPlotParameters(object=this, fcn="plotXY", what=what, slide=slide, 
                                                      xlog=xlog, ylog=ylog);
  if (is.null(cmd)) {  
    if (length(pch) == 1 && pch == "auto")
      pch <- 176;
    plot(x, y, cex=cex, pch=pch, col=col, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...);

    Device$putTimestamp();
    Device$putDataset();
    Device$putAuthor();
    putSlide(this);
  } else if (cmd == "points") {
    if (length(pch) == 1 && pch == "auto")
      pch <- 176;
    points(x, y, cex=cex, pch=pch, col=col, ...);
  } else if (cmd == "highlight") {
    if (length(pch) == 1 && pch == "auto")
      pch <- 20;
    points(x, y, cex=cex, pch=pch, col=col, ...);
  } else if (cmd == "text") {
    text(x, y, cex=cex, col=col, ...);
  } else if (cmd == "lowessline") {
    # Line types:
    # (0=blank, 1=solid, 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash)
    # Plot the lowess curve for ALL pins (library: sma)
    idx <- !(is.na(x) | is.na(y) | is.infinite(x) | is.infinite(y));

    line <- approx(lowess(x[idx], y[idx], f=f), ties=mean);
    returnValue <- line;

    lines(line, col=col, ...);

    if (gridwise) {
      # Plot the lowess curve for EACH pin individually
  
      x.all[!include] <- NA;
      y.all[!include] <- NA;
  
      layout <- getLayout(this);
      if (is.null(layout))
  	throw("Layout is not specified. Don't know how to plot gridwise.");
  
      if (is.null(col)) col <- rainbow(layout$ngrid.c);
      col <- rep(col, length.out=nbrOfGrids(layout));
  
      # Spot indices where each pin starts.
      # Comment: This assumes that the gene expression values are listed
      #          gridwise: 1, 2, 3, ..., #grids.
      nrow <- layout$ngrid.r;
      ncol <- layout$ngrid.c;
      npin <- nrow*ncol;
      pin <- c(0, rep(layout$nspot.r*layout$nspot.c, npin) * (1:npin));
      pin.lty <- matrix(1:nrow, nrow=nrow, ncol=ncol, byrow=TRUE);
      pin.lty <- rep(pin.lty, length.out=npin);

      idx <- !(is.na(x.all) | is.infinite(x.all) | 
               is.na(y.all) | is.infinite(y.all)  );
      nbrOfSpots <- 1:nbrOfSpots(this);
      line <- list()
      for (j in 1:npin) {
  	    pin.idx <- ((pin[j] + 1) <= nbrOfSpots) & (nbrOfSpots <= pin[j+1]);
        index <- idx & pin.idx;
  	    x <- x.all[index];
  	    y <- y.all[index];
        line[[j]] <- approx(lowess(x, y, f=f), ties=mean);
      }
      for (j in 1:npin)
        lines(line[[j]], col="gray", lty=1, ...);
      for (j in 1:npin)
        lines(line[[j]], col=col[j], lty=pin.lty[j], ...);
  	
      layout$put(x="99%", y="1%", col="auto", lty="auto");
    }
  } else {
    throw("Unkown plot type: ", cmd);
  }

  invisible(returnValue);
})



#########################################################################/**
# @RdocMethod plotSpatial
#
# @title "Plots a spatial representation of one field"
#
# \description{
#  @get "title".
#  It is recommended to use the \code{plot} function instead of calling
#  this method explicitly (see @seemethod "plot").
# }
#
# @synopsis
#
# \arguments{
#  \item{what}{What to plot. Any field that can be retrieved by \code{extract},
#   is accepted.}
#  \item{slide}{The slide to be plotted.}
#  \item{include}{The indices of the spots that should be included. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are considered.
#   If @NULL all spots are considered.}
#  \item{exclude}{The indices of the spots that should be excluded. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are excluded. If @NULL no spots are excluded.}
#  \item{col}{The color(s) to be used for the plotted spots, i.e. for the
#   spots \emph{after} inclusion and exclusion. If the value is
#   \code{"redgreen"} a red to green palette is used.}
#  \item{...}{Common arguments accepted by most plot functions.
#   For more information see @see "graphics::par" and @see "graphics::plot".}
#  \item{cex}{For internal use only! See above.}
#  \item{pch}{For internal use only! See above.}
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   subplots(4)
#   plotSpatial(ma)                   # Spatial plot of log ratios before.
#   normalizeWithinSlide(ma, "p")     # Printtipwise lowess normalization.
#   plotSpatial(ma)                   # Spatial plot of log ratios after.
#
#   plotSpatial(ma, include=(abs(ma$M) > 2))
#   points(ma, include=(abs(ma$M) > 2), col="red")
# }
#
# @author
#
# \seealso{
#   @seemethod "plot".
#   @seemethod "plotXY".
#   @seemethod "highlight".
#   @seemethod "text".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("plotSpatial", "MicroarrayData", function(this, what, slide=1, include=NULL, exclude=NULL, col="auto", xlab=NULL, ylab="", axes=FALSE, xaxs="i", yaxs="i", pch="auto", grid=TRUE, log=NULL, ..., cex=NULL, style=NULL) {
  # To be removed in August 2005. /HB 2005-05-04
  if (!exists("image270", mode="function"))
    image270 <- image.matrix;

  cmd <- NULL;
  if (!is.null(style) && is.element(style, c("points", "highlight", "text"))) {
    cmd <- style;
    lastPlot <- Device$getPlotParameters();
    what <- lastPlot$what;
    slide <- lastPlot$slide;
  }

  slide <- validateArgumentSlide(this, slide=slide);

  colWhat <- what;
  if (identical(what, c("M", "A")))
    what <- c("A","M");

  if (identical(what, c("A", "M")))
    what <- "M";

  if (length(what) != 1 || !is.character(what))
    throw("Argument 'what' must be a single field name.");

  X <- what;

  setView(this, MicroarrayArray$SPOT.SLIDE.VIEW);

  include <- which(getInclude(this, include, exclude, slide=slide));
  
  if (length(col) == 1) {
    color <- col;
    col <- rep(NA, length.out=nbrOfSpots(this));
    if (!is.numeric(color) && substring(color,1,1) != "#" && !is.element(color, colors())) {
      # Generate the colors automatically using the value of specified field.
      col[include] <- getColors(this, what=colWhat, slide=slide, include=include, palette=color, log=log);
    } else {
      col[include] <- color;
    }
  }

  layout <- getLayout(this);

  Device$setPlotParameters(object=this, fcn="plotSpatial", what=what,
                           slide=slide);

  if (is.null(cmd)) {
    if (is.null(xlab)) xlab <- getLabel(this, X);
    if (is.expression(xlab)) xlab <- xlab[[1]];
    colorMap <- sort(unique(col));
    col <- match(col, colorMap);
    col <- rep(col, length.out=nbrOfSpots(layout))
    col <- toXYMatrix(layout, col)
    x <- 1:(ncol(col)+0);
    y <- 1:(nrow(col)+0);
    image270(x=x, y=y, z=col, col=colorMap, xaxs=xaxs, yaxs=yaxs, xlab=xlab, ylab=ylab, bty="n", axes=axes, ...);

    Device$putTimestamp();
    Device$putDataset();
    Device$putAuthor();
    putSlide(this);
  } else if (cmd == "points") {
    colorMap <- sort(unique(col));
    col <- match(col, colorMap);
    col <- rep(col, length.out=nbrOfSpots(layout))
    col <- toXYMatrix(layout, col)
    x <- 1:(ncol(col)+0);
    y <- 1:(nrow(col)+0);
    image270(x=x, y=y, z=col, col=colorMap, xaxs=xaxs, yaxs=yaxs, xlab=xlab, ylab=ylab, bty="n", axes=axes, add=TRUE, ...);
  } else if (cmd == "highlight") {
    if (length(pch) == 1 && pch == "auto") pch <- 20
    col <- col[include];
    positions <- getPosition(layout, include);
    xx <- positions[,2];
    yy <- positions[,1];
    points(xx, nbrOfRows(layout)+1-yy, cex=cex, col=col, pch=pch, ...);
  } else if (cmd == "text") {
    col <- col[include]
    positions <- getPosition(layout, include);
    xx <- positions[,2];
    yy <- positions[,1];
    text(xx, nbrOfRows(layout)+1-yy, cex=cex, col=col, ...);
  } else {
    throw("Unknown plot type: ", cmd);
  }

  if (grid) {
    # Draw the grids and the border around the plot
    gr <- layout$ngrid.r;
    gc <- layout$ngrid.c;
    sr <- layout$nspot.r;
    sc <- layout$nspot.c;
    box();
    abline(h = ((gr-1):1)*sr + 0.5);
    abline(v = (1:(gc-1))*sc + 0.5);
  }
})




#########################################################################/**
# @RdocMethod plotDensity
#
# @title "Plots the empirical density distribution of a field"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{what}{What to plot. Any field that can be retrieved by
#   \code{extract} is accepted.}
#  \item{slides}{The slides to be plotted. Each slide will be plotted
#   as a separate curve.}
#  \item{log}{The log base to be used for taking the logarithm of the
#   data before generating the plot. If @NULL the non-logarithm data
#   is plotted.}
#  \item{xlim,ylim}{@character @vector of length 2. The x and y limits.}
#  \item{xlab,ylab}{@character string for labels on x and y axis.}
#  \item{col}{The color(s) of the curves.}
#  \item{add}{If @TRUE, the curves are plotted in the current plot,
#   otherwise a new is created.}
#  \item{legend}{If @TRUE, a legend is added, otherwise not.}
#  \item{...}{Common arguments accepted by most plot functions.
#   For more information see @see "graphics::par" and @see "graphics::plot".}
# }
#
# @author
#
# \seealso{
#   @seemethod "plot".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("plotDensity", "MicroarrayData", function(this, what, slides=NULL, log=NULL, xlim=NULL, ylim=NULL, xlab=what, ylab="density (integrates to one)", col=NULL, add=FALSE, legend=!add, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(col))
    col <- slides;

  # Roll out all vectors
  col <- rep(col, length.out=max(slides));

  if (add & is.null(xlim))
    xlim <- par("usr")[1:2];

  # Generate all densities first and figure out the plot limits.
  ds <- list();
  xlimDef <- c(NA,NA);
  ylimDef <- c(0,NA);
  for(kk in slides) {
    x <- this[[what]][,kk];
    x <- na.omit(x);
    if (!is.null(log)) {
      x <- log(x, base=log);
      x <- x[is.finite(x)];
    }
    # If xlim is specified, cut the data first, do get a higher resolution of 
    # the density in that region to be shown.
    if (!is.null(xlim))
      x <- x[x >= xlim[1] & x <= xlim[2]];
    d <- density(x);
    rm(x);
    ds[[kk]] <- d;
    xlimDef <- range(c(xlimDef, range(d$x, na.rm=TRUE)), na.rm=TRUE);
    ylimDef <- range(c(ylimDef, range(d$y, na.rm=TRUE)), na.rm=TRUE);
    rm(d);
  }
  gc();

  if (is.null(xlim))
    xlim <- xlimDef;
  if (is.null(ylim))
    ylim <- ylimDef;

  if (add == FALSE)
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...);
  for(kk in slides) {
    lines(ds[[kk]], col=col[kk], ...);
  }

  if (identical(legend, TRUE)) {
    # Put a legend in the upper right corner
    usr <- par("usr")
    x <- usr[2]
    y <- usr[4]
    pmt <- getSlideNames(this);
    if (is.null(pmt))
      pmt <- slides;
    legend(x=x,y=y, legend=pmt, fill=col, xjust=1, yjust=1, cex=0.7);
  }
}) # plotDensity()





#########################################################################/**
# @RdocMethod plotPrintorder
#
# @title "Plots the data as a time series in the order it was printed"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{field}{The field to be plotted. Any field that can be retrieved
#       by \code{extract}, is accepted.}
#  \item{slides}{The slide(s) to be plotted. If @NULL, all slides are
#       plotted.}
#  \item{include}{The indices of the spots that should be included. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are considered.
#   If @NULL all spots are considered.}
#  \item{exclude}{The indices of the spots that should be excluded. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are excluded. If @NULL no spots are excluded.}
#  \item{col}{The color(s) to be used for the plotted spots, i.e. for the
#   spots \emph{after} inclusion and exclusion. If the value is
#   \code{"redgreen"} a red to green palette is used.}
#  \item{...}{Common arguments accepted by most plot functions.
#   For more information see @see "graphics::par" and @see "graphics::plot".}
# }
#
# \details{
#   If a vector of colors (\code{col}), symbol sizes (\code{cex}) or
#   symbol types (\code{pch}) are specified, such vectors will be reordered
#   according to the print order \emph{before} being applied.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   plotPrintorder(ma, "M")       # Printorder plot of log ratios before.
#   normalizeWithinSlide(ma, "p") # Printtipwise lowess normalization.
#   plotPrintorder(ma, "M")       # Printorder plot of log ratios after.
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
#   @seemethod "plot".
#   @seemethod "plotXY".
#   @seemethod "highlight".
#   @seemethod "text".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("plotPrintorder", "MicroarrayData", function(this, what, slides=NULL, include=NULL, exclude=NULL, cex="auto", col="auto", pch="auto", xlab="in print order", ylab=NULL, ylog=NULL, groupfcn=NULL, breakpoints=NULL, ..., f=200/nbrOfSpots(this), style=NULL) {

  returnValue <- NULL;
  cmd <- NULL;
  if (!is.null(style) && is.element(style, c("points", "lowessline", "highlight", "text"))) {
    cmd <- style;
    lastPlot <- Device$getPlotParameters();
    what <- lastPlot$what;
    slides <- lastPlot$slides;
    ylog <- lastPlot$ylog;
  }

  slides <- validateArgumentSlides(this, slides=slides);

  if (length(what) != 1 || !is.character(what))
    throw("Argument 'what' must be a single field name.");

  if (!hasLayout(this))
    throw("Can not plot the data in print order since the layout is not specified.");

  layout <- getLayout(this);
  printorderMatrix <- toPrintorderMatrix(layout);
  printorder       <- as.vector(printorderMatrix);
  
  setView(this, MicroarrayArray$SPOT.SLIDE.VIEW);

  # Get which spots to include
  include <- which(getInclude(this, include, exclude, slides=slides));

  # Match these spots towards the printorder indices.
  match <- na.omit(match(include, printorder));

  # Set all non-includes to NA
  if (length(match) > 0)
    printorder[-match] <- NA;

  xs <- c(); ys <- c();
  for (slide in slides) {
    ystmp <- getSpotSlideValues(this[[what]], spots=printorder, slides=slide);

    # Now ys is a vector, but if we would to e.g. take the mean of all groups
    # We have to rearrange ys back to a matrix where each column represents a group.
    if (!is.null(groupfcn)) {
      ysMatrix <- matrix(ystmp, nrow=nrow(printorderMatrix));
      # Apply the wanted function on each column, i.e. ys will become shorter
      ystmp <- apply(ysMatrix, MARGIN=2, FUN=groupfcn);
    }
    ys <- c(ys, ystmp);
    xs <- c(xs, 1:length(ystmp));
  }
  
  if (!is.null(ylog) && ylog > 0)
    ys <- log(ys, base=ylog);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # C o l o r   s e t u p
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!is.null(col)) {
    if (identical(col, "auto")) {
      palette <- col;
      if (isFieldColorable(this, what)) {
  	col <- c();
  	for (slide in slides) {
  	  col0 <- getColors(this, what=what, slide=slide, palette=palette);
  	  col <- cbind(col, col0);
  	}
      } else {
  	col <- NULL;
      }
    } else {
      col <- as.matrix(col);
    }
  }
  if (length(col) > 1)
    col <- col[printorder,];

  if (!is.null(cex)) {
    if (cex == "auto") cex <- NULL; 
  } else if (length(cex) > 1) {
    cex <- cex[include];
    cex <- as.matrix(cex);
    cex <- cex[printorder,];
  }

  if (!is.null(pch)) {
    if (length(pch) > 1) {
      pch <- pch[include];
      pch <- as.matrix(pch);
      pch <- pch[printorder,];
    }
  }

  Device$setPlotParameters(object=this, fcn="plotPrintorder", what=what,
                                                slides=slides, ylog=ylog);

  if (is.null(cmd)) {
    if (length(pch) == 1 && pch == "auto")
      pch <- 176;

    if (is.null(ylab)) ylab <- getLabel(this, what);
    if (is.expression(ylab)) ylab <- ylab[[1]];
  
    plot(xs, ys, cex=cex, col=col, pch=pch, xlab=xlab, ylab=ylab, ...);
    if (!is.null(breakpoints)) {
      abline(v=breakpoints, col="gray")
    }

    Device$putTimestamp();
    Device$putDataset();
    Device$putAuthor();
    putSlide(this);
  } else if (cmd == "points") {
    if (length(pch) == 1 && pch == "auto")
      pch <- 176;
    points(xs, ys, cex=cex, col=col, pch=pch, ...);
  } else if (cmd == "highlight") {
    if (length(pch) == 1 && pch == "auto")
      pch <- 20;
    points(xs, ys, cex=cex, col=col, pch=pch, ...);
  } else if (cmd == "text") {
    text(xs, ys, cex=cex, col=col, pch=pch, ...);
  } else if (cmd == "lowessline") {
    idx <- !(is.na(xs) | is.na(ys) | is.infinite(xs) | is.infinite(ys));
    line <- lowess(xs[idx], ys[idx], f=f);
#    line <- approx(line, ties=mean);
    returnValue <- line;
    lines(line, col="white", lwd=2);
    lines(line, col=col, ...);
  } else {
    throw("Unknown plot type: ", cmd);
  }

  if (!is.null(returnValue)) invisible(returnValue);
})





setMethodS3("plotDiporder", "MicroarrayData", function(this, what, slides=NULL, include=NULL, exclude=NULL, cex="auto", col="auto", pch="auto", xlab=paste("in dip order (", nbrOfGrids(getLayout(this)), " spots/dip)", sep=""), ylab=NULL, ylog=NULL, groupfcn=NULL, breakpoints=NULL, ..., f=200/nbrOfSpots(this), style=NULL) {

  returnValue <- NULL;
  cmd <- NULL;
  if (!is.null(style) && is.element(style, c("points", "lowessline", "highlight", "text"))) {
    cmd <- style;
    lastPlot <- Device$getPlotParameters();
    what <- lastPlot$what;
    slides <- lastPlot$slides;
    ylog <- lastPlot$ylog;
  }

  if (length(what) != 1 || !is.character(what))
    throw("Argument 'what' must be a single field name.");

  slides <- validateArgumentSlides(this, slides=slides);

  if (!hasLayout(this))
    throw("Can not plot the data in print order since the layout is not specified.");

  setView(this, MicroarrayArray$SPOT.SLIDE.VIEW);

  # Get which spots to include
  include <- which(getInclude(this, include, exclude, slides=slides));

  layout <- getLayout(this);
  dips <- getPrintdipGroups(layout);
  spots <- getSpots(dips);

  # The data to be plotted
  data <- getSpotSlideValues(this[[what]], slides=slides);

  # Keep only data to be included, all others should be NA's
  data[-include] <- NA;
  
  if (is.null(groupfcn)) groupfcn <- median;

  xs <- 1:length(spots);
  ys <- rep(NA, length(spots));
  for (k in xs) {
    x <- data[spots[[k]]]; x <- x[!is.na(x)]; x <- x[!is.infinite(x)];
    if (length(x) > 0)
      ys[k] <- groupfcn(x);
  }

  if (!is.null(ylog) && ylog > 0)
    ys <- log(ys, base=ylog);

  if (!is.null(col) && col == "auto") {
    palette <- col;
    if (what %in% c("R", "Rb")) {
      if (max(ys, na.rm=TRUE) > 100)
        col <- Colors$getRed(log(ys,2), x.range=c(0,16))
      else
        col <- Colors$getRed(ys, x.range=c(0,16));
    } else if (what %in% c("G", "Gb")) {
      if (max(ys, na.rm=TRUE) > 100)
        col <- Colors$getGreen(log(ys,2), x.range=c(0,16))
      else
        col <- Colors$getGreen(ys, x.range=c(0,16));
    } else if (what == "M") {
      A <- rep(10, length(ys));
      col <- MicroarrayData$createColors(A, ys);
    } else if (what == "A") {
      col <- Colors$getGray(ys, x.range=c(0,16));
    } else {
      col <- Colors$getGray(ys, x.range=c(0,16));
    }
  }

  if (!is.null(cex) && cex == "auto") cex <- NULL; 

  if (length(cex) > 1) cex <- cex[include];
  if (length(pch) > 1) pch <- pch[include];


  Device$setPlotParameters(object=this, fcn="plotDiporder", what=what,
                                              slides=slides, ylog=ylog);

  if (is.null(cmd)) {
    if (length(pch) == 1 && pch == "auto")
      pch <- 176;

    if (is.null(ylab)) ylab <- getLabel(this, what);
    if (is.expression(ylab)) ylab <- ylab[[1]];
  
    plot(xs, ys, cex=cex, col=col, pch=pch, xlab=xlab, ylab=ylab, ...);
    if (!is.null(breakpoints)) {
      abline(v=breakpoints, col="gray")
    }
    Device$putTimestamp();
    Device$putDataset();
    Device$putAuthor();
    putSlide(this);
  } else if (cmd == "points") {
    if (length(pch) == 1 && pch == "auto")
      pch <- 176;
    points(xs, ys, cex=cex, col=col, pch=pch, ...);
  } else if (cmd == "highlight") {
    if (length(pch) == 1 && pch == "auto")
      pch <- 20;
    points(xs, ys, cex=cex, col=col, pch=pch, ...);
  } else if (cmd == "text") {
    text(xs, ys, cex=cex, col=col, pch=pch, ...);
  } else if (cmd == "lowessline") {
    idx <- !(is.na(xs) | is.na(ys) | is.infinite(xs) | is.infinite(ys));
    line <- lowess(xs[idx], ys[idx], f=f);
#    line <- approx(line, ties=mean);
    returnValue <- line;
    lines(line, col="white", lwd=2);
    lines(line, col=col, ...);
  } else {
    throw("Unknown plot type: ", cmd);
  }

  if (!is.null(returnValue)) invisible(returnValue);
}, private=TRUE)







setMethodS3("plotReplicates", "MicroarrayData", function(this, what, gene=1, slides=NULL, pch=176, hline=median, ylog=NULL, new=TRUE, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  returnValue <- NULL;
 
  layout <- getLayout(this);
  if (hasReplicates(layout)) {
    reps <- getReplicates(layout);
    index <- getSpot(reps, gene);
  } else
    index <- gene;

  field <- this[[what]];
  ys <- field[index,slides];

  if (!is.null(ylog) && ylog > 0)
    ys <- log(ys, base=ylog);

  Device$setPlotParameters(object=this, fcn="plotReplicates", what=what, gene=gene, ...);

  if (new == TRUE) {
    plot(seq(ys), ys, pch=pch, ...)
    Device$putTimestamp();
    Device$putDataset();
    Device$putAuthor();
    putGene(this);
  } else {
    points(seq(ys), ys, pch=pch, ...);
  }

  if (is.matrix(ys) && ncol(ys) > 1) {
    # Has within-slide replicates

    # Draw lines between within-slide replicates
		for (k in seq(ncol(ys))) {
 			xs <- seq(2*k-1, length=length(index));
 			lines(xs, ys[xs], ...);
    }

  	# Draw lines between between-slide replicates
 		for (k in seq(from=1, to=ncol(ys)-1)) {
 			xs <- seq(2*k, length=length(index));
 			lines(xs, ys[xs], lty=2, ...);
  	}
  } else {
    # No within-slide replicates
  	# Draw lines between between-slide replicates
  	if (length(ys) > 1) {
 			xs <- seq(1, length=length(ys));
 			lines(xs, ys[xs], lty=2, ...);
  	}
  }

  if (!is.null(hline)) 
    abline(h=hline(na.omit(ys)), lty=3, ...);

  if (!is.null(returnValue)) invisible(returnValue);
})





#########################################################################/**
# @RdocMethod boxplot
#
# @title "Plots a boxplot"
#
# \description{
#  Creates a box plot for any data over all grids or over all slides. If
#  the argument \code{slide} is given, an boxplot over the print-tips will
#  be generated, others an across-slides boxplot will be generated.
# }
#
# @synopsis
#
# \arguments{
#  \item{what}{What data variable to do a boxplot for.}
#  \item{slides}{The slide(s) that should be used to generate this plot.
#     If @NULL data from all slides are considered.}
#   \item{include}{The indices of the spots that should be included. 
#    If it is instead a name of one or more flags, the spots which have been
#    flagged with these flags are considered.
#    If @NULL all spots are considered.}
#   \item{exclude}{The indices of the spots that should be excluded. 
#    If it is instead a name of one or more flags, the spots which have been
#    flagged with these flags are excluded.}
#  \item{groupBy}{The data can be grouped by either \code{"slide"}
#    or \code{"printtip"}.}
#  \item{gridwise}{(Deprecated) 
#     If @TRUE a within-slide (across-grid) boxplot 
#     is generated, otherwise an across-slide boxplot is generated.}
#  \item{names}{A string vector of names for each grid/slide.
#   If @NULL the grids will be named 1,2,3, etc.}
#  \item{xlab, ylab}{Label on the x-axis (y-axis). If @NULL 
#   the label will be automatically set.}
#  \item{las}{Rotation style of the axis labels. Default value is 3 
#   (always vertical). For more information see @see "graphics::par".}
#  \item{...}{Common arguments accepted by most plot functions.
#   For more information see @see "graphics::par" and @see "graphics::plot".}
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   boxplot(ma, groupBy="printtip") # Log-ratios for printtips on all slides
#   boxplot(ma)                     # Log-ratios for all slides.
#   boxplot(ma, slides=c(1,3,5))    # Log-ratios for slides 1, 3 and 5.
# }
#
# @author
#
# \seealso{
#   @see "graphics::boxplot".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("boxplot", "MicroarrayData", function(x, what, style=NULL, 
                   slides=NULL, include=NULL, exclude=NULL, groupBy="slide", 
                   gridwise=FALSE, 
                   names=NULL, xlab=NULL, ylab=NULL, log=NULL, las=3, 
                   col=NULL, ...) {
  # To please R CMD check...
  this <- x;

  cmd <- NULL;
  if (!is.null(style) && is.element(style, c("points", "highlight", "text"))) {
    cmd <- style;
    lastPlot <- Device$getPlotParameters();
    what <- lastPlot$what;
    slides <- lastPlot$slides;
    groupBy <- lastPlot$groupBy;
  }

  if (length(what) != 1)
    throw("Argument 'what' in call to boxplot must be a single field name.");

  slideMissing <- is.null(slides);
  if (is.null(slides))
    slides <- seq(from=1, to=ncol(this[[what]]));
  nbrOfSlides <- length(slides);

  data <- extract(this, field=what, slide=slides);
  data <- matrix(as.matrix(data), ncol=nbrOfSlides);

  include <- getInclude(this, include=include, exclude=exclude, slides=slides);
  include <- which(include);

#  data <- data[include]; # <-- Should this be data[include,] instead?

  # Ad hoc solution below and this should be fixed. /HB 2001-09-30
  data <- matrix(as.matrix(data), ncol=nbrOfSlides);
  
  if (is.null(ylab))
    ylab <- this$getLabel(what);
  if (is.expression(ylab))
    ylab <- ylab[[1]];

  if (!is.null(log)) {
    data <- log(data, base=log);
    ylab <- substitute(log[YY](XX), list(XX=ylab, YY=log));
  }

  data[is.infinite(data)] <- NA;

  if (gridwise)
    groupBy <- "printtips";
  if (is.null(groupBy))
    groupBy <- "slides";

  setXlab <- (is.null(xlab));

  layout <- getLayout(this);
  groups <- NULL;
  if (inherits(groupBy, "LayoutGroups")) {
    groups <- groupBy;
    if (setXlab) {
      s <- data.class(groupBy);
      s <- gsub("Groups$", "", s);
      s <- unlist(strsplit(s, ""));
      upper <- which(s == toupper(s));
      upper <- upper[upper > 0];
      xlab <- c();
      for (k in seq(along=s)) {
        if (s[k] == toupper(s[k]) && k > 1) {
          xlab <- c(xlab, " ");
          s[k] <- tolower(s[k]);
        }
        xlab <- c(xlab, s[k]);
      }
      xlab <- paste(xlab, collapse="");
    }
  } else if (is.character(groupBy)) {
    groupBy <- tolower(groupBy);
    if (groupBy == "plate") {
      groups <- getPlateGroups(layout);
      if (setXlab)
  	xlab <- "Plate";
    } else if (groupBy == "printdip") {
      groups <- getPrintdipGroups(layout);
      if (setXlab)
  	xlab <- "Print dip";
    } else if (groupBy == "printtip") {
      groups <- getPrinttipGroups(layout);
      if (setXlab)
  	xlab <- "Print tip";
    } else if (groupBy == "sliderow") {
      groups <- getSlideRowGroups(layout);
      if (setXlab)
  	xlab <- "Slide row";
    } else if (groupBy == "slidecolumn") {
      groups <- getSlideColumnGroups(layout);
      if (setXlab)
  	xlab <- "Slide column";
    } else if (groupBy == "slide") {
      groups <- NA;
      if (setXlab)
  	xlab <- "Slide";
    }
  }

  if (is.null(groups))
      throw("Unknown value of argument 'groupBy': ", groupBy);

  if (setXlab) {
    if (nbrOfSlides(this) > 1) {
      if (slideMissing)
    	xlab <- paste(xlab, "- all slides")
      else
    	xlab <- paste(xlab, "- slides", paste(slides, collapse=", "));
    }
  }

  if (inherits(groups, "LayoutGroups")) {
    spots <- getSpots(groups);
    if (is.null(names))
      names <- getNames(groups);
    rm(groups);
    n <- unlist(lapply(spots, FUN=length));
    df <- matrix(NA, nrow=nbrOfSlides*max(n), ncol=length(spots));
    for (k in seq(along=spots)) {
      idx <- intersect(spots[[k]], include);
      if (length(idx) > 0) {
        value <- as.vector(data[idx,]);
        df[1:length(value),k] <- value;
      }
    }
    data <- df;
    rm(df);
  } else {
    if (is.null(names) && !is.null(slides))
      names <- slides;
  }

  if (is.null(col)) {
    avg <- apply(data, MARGIN=2, FUN=median, na.rm=TRUE);
    if (what == "M") {
      tmp <- clone(this);
      avg <- as.matrix(avg);
      tmp[["A"]] <- matrix(13, nrow=nrow(avg), ncol=ncol(avg));
      tmp[[what]] <- avg;
      col <- getColors(tmp, what=what);
    } else if (what == "A") {
      col <- Colors$getGray(avg);
    } else {
    }
  }

  data <- as.data.frame(data);

  Device$setPlotParameters(object=this, fcn="boxplot", what=what,
    groupBy=groupBy, slides=slides, groupBy=groupBy, log=log);

  if (is.null(cmd) || identical(cmd, "boxplot")) {
    boxplot(data, names=names, las=las, col=col, ...); 
    title(xlab=xlab, ylab=ylab);

    Device$putTimestamp();
    Device$putDataset();
    Device$putAuthor();
    putSlide(this);
  } else if (identical(cmd, "highlight")) {
    x <- seq(ncol(data));
    x <- matrix(x, nrow=nrow(data), ncol=ncol(data), byrow=TRUE);
    points(x, as.matrix(data), col=col, ...);
  }

  invisible(data);
}) # boxplot()



setMethodS3("qqnorm", "MicroarrayData", function(y, what, slides=NULL, include=NULL, exclude=NULL, pch="auto", main="Normal Q-Q Plot", xlab="Quantiles of standard normal", ylab=what, ..., style=NULL) {
  # To please R CMD check...
  this <- y;

  cmd <- NULL;
  if (!is.null(style) && is.element(style, c("points", "highlight", "text"))) {
    cmd      <- style;
    lastPlot <- Device$getPlotParameters();
    what     <- lastPlot$what;
    slides   <- lastPlot$slides;
    qq       <- lastPlot$qq;
  } else {
    qq <- NULL;
  }

  if (length(what) != 1)
    throw("Argument 'what' in call to boxplot must be a single field name: ", what);

  slides <- validateArgumentSlides(this, slides=slides);

  include <- getInclude(this, include=include, exclude=exclude, slides=slides);

  if (is.null(qq)) {
    data <- extract(this, field=what, slide=slides);
    nbrOfSlides <- if (is.null(slides)) nbrOfSlides(this) else length(slides);
    # Ad hoc solution below and this should be fixed. /HB 2001-09-30
    data <- matrix(as.matrix(data), ncol=nbrOfSlides);
    data <- data[include,];
    qq <- qqnorm(data, plot.it=FALSE);
    x <- qq$x;
    y <- qq$y;
  } else {
    x <- qq$x[include];
    y <- qq$y[include];
  }

  Device$setPlotParameters(object=this, fcn="qqnorm", what=what,
    slides=slides, qq=qq);

  if (is.null(cmd)) {  
    if (length(pch) == 1 && pch == "auto")
      pch <- 176
    plot(x,y, pch=pch, xlab=xlab, ylab=ylab, ...);
    title(main=main);
    Device$putTimestamp();
    Device$putDataset();
    Device$putAuthor();
    putSlide(this);
  } else if (cmd == "points") {
    if (length(pch) == 1 && pch == "auto")
      pch <- 176
    points(x, y, pch=pch, ...);
  } else if (cmd == "highlight") {
    if (length(pch) == 1 && pch == "auto")
      pch <- 20
    points(x, y, pch=pch, ...);
  } else if (cmd == "text") {
    text(x, y, ...);
  }

  invisible(qq);
}) # qqnorm()



#########################################################################/**
# @RdocMethod hist
#
# @title "Plots a histogram"
#
# \description{
#  Creates a histogram plot for any data variable over all grids or over
#  all slides. If the argument \code{slide} is given, an boxplot over the print-tips will
#  be generated, others an across-slides boxplot will be generated.
# }
#
# @synopsis
#
# \arguments{
#   \item{what}{What data variable to do a histogram for.}
#   \item{slides}{The slide(s) that should be used to generate this plot.
#    If @NULL data from all slides are considered.}
#   \item{include}{The indices of the spots that should be included. 
#    If it is instead a name of one or more flags, the spots which have been
#    flagged with these flags are considered.
#    If @NULL all spots are considered.}
#   \item{exclude}{The indices of the spots that should be excluded. 
#    If it is instead a name of one or more flags, the spots which have been
#    flagged with these flags are excluded. If @NULL no spots are excluded.}
#   \item{xlab, ylab}{Label on the x-axis (y-axis). If @NULL 
#    the label will be automatically set.}
#   \item{log}{The log base to be used for taking the logarithm of the data 
#     before generating the histogram. If @NULL the non-logarithm data
#     is plotted.}
#   \item{...}{Common arguments accepted by most plot functions.
#     For more information see @see "graphics::par" and @see "graphics::plot".}
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   hist(ma, what="M", slide=1)
# }
#
# @author
#
# \seealso{
#   @see "graphics::boxplot".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("hist", "MicroarrayData", function(x, what, slides=NULL, include=NULL, exclude=NULL, xlab=NULL, log=NULL, ...) {
  # To please R CMD check...
  this <- x;

  slides <- validateArgumentSlides(this, slides=slides);

  if (length(what) != 1)
    throw("Argument 'what' in call to hist must be a single field name.");

  data <- extract(this, field=what, slide=slides);
  nbrOfSlides <- if (is.null(slides)) nbrOfSlides(this) else length(slides);
  data <- matrix(as.matrix(data), ncol=nbrOfSlides);
 
  if (is.null(xlab))
    xlab <- this$getLabel(what);
  if (is.expression(xlab))
    xlab <- xlab[[1]];

  include <- getInclude(this, include=include, exclude=exclude, slides=slides);
  data <- as.vector(data[include]);
  if (!is.null(log))
    data <- log(data, base=log);
  
  # Exclude all NA/NaN/Inf's
  ok <- !(is.na(data) | is.infinite(data));
  data <- data[ok];

  Device$setPlotParameters(object=this, fcn="hist", what=what,
    slides=slides, log=log);

  hist(data, xlab=xlab, ...);
  Device$putTimestamp();
  Device$putDataset();
  Device$putAuthor();
  putSlide(this);
})



############################################################################
# Plot annotations
############################################################################

#########################################################################/**
# @RdocMethod highlight
#
# @title "Highlight spots in last plot"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{include}{The indices of the spots that should be included. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are considered.
#   If @NULL all spots are considered.}
#  \item{exclude}{The indices of the spots that should be excluded. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are excluded.
#   If @NULL no spots are excluded.}
#  \item{pch}{The point style to be used. Default value is 176 (a small open circle). For more information see @see "graphics::par".}
#  \item{cex}{The scale factor to be used. If @NULL the system default scale factor will be used. For more information see @see "graphics::par".}
#  \item{col}{The color(s) to be used for the plotted spots, i.e. for the
#   spots \emph{after} inclusion and exclusion.}
#  \item{...}{Common arguments accepted by most plot functions.
#   For more information see @see "graphics::par" and @see "graphics::plot".}
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   # Get all spots in the grid at row 1 and at column 2.
#   layout <- getLayout(ma)
#   idx <- getIndices(layout, 1,2, NULL,NULL)
#
#   # Plot the data both highlighted and non-highlighted.
#   layout(matrix(1:4, ncol=2, byrow=TRUE))
#   plot(ma)
#   plotSpatial(ma)
#   plot(ma)
#   highlight(ma, idx, col="purple")  # Highlight all spots in grid (1,2).
#   plotSpatial(ma)
#   highlight(ma, idx, col="purple")  # Highlight all spots in grid (1,2).
# }
#
# @author
#
# \seealso{
#   @seemethod "text".
#   @see "graphics::points".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("highlight", "MicroarrayData", function(this, include=NULL, exclude=NULL, col="black", ...) {
  lastPlot <- Device$getPlotParameters();
  if (length(lastPlot) == 0)
    throw("To call points() of a MicroarrayData class one must first call one of the plot functions of a MicroarrayData class that creates a new plot must be called first.");
  fcn <- get(lastPlot$fcn, mode="function");
  fcn(this, style="highlight", include=include, exclude=exclude, col=col, ...);
})



setMethodS3("points", "MicroarrayData", function(x, what=NULL, include=NULL, exclude=NULL, col="auto", ...) {
  # To please R CMD check...
  this <- x;

  lastPlot <- Device$getPlotParameters();
  if (length(lastPlot) == 0)
    throw("To call points() of a MicroarrayData class one must first call one of the plot functions of a MicroarrayData class that creates a new plot must be called first.");
  fcn <- get(lastPlot$fcn, mode="function")
  if (!is.null(what)) {
    searches <- c("vs");
    match <- which(apply(as.matrix(searches), MARGIN=1, 
                        FUN=function(x) regexpr(x, what)) != -1);
    what <- as.character(what);
    what <- rev(unlist(strsplit(what, searches[match])));
    what <- what[nchar(what) != 0];
  }
  fcn(this, style="points", what=what, include=include, exclude=exclude, col=col, ...);
})


#########################################################################/**
# @RdocMethod text
#
# @title "Puts labels to the spots in last plot"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{include}{The indices of the spots that should be included. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are considered.
#   If @NULL all spots are considered.}
#  \item{exclude}{The indices of the spots that should be excluded. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are excluded.
#   If @NULL no spots are excluded.}
#  \item{labels}{Labels of \emph{all} the data points. If @NULL first the labels of the belonging \link{Layout} object is retrieved. If there is no Layout object, the data points are labeled with indices of the plotted spots.}
#  \item{col}{The color(s) to be used for the plotted spots, i.e. for the
#   spots \emph{after} inclusion and exclusion.}
#  \item{pos}{A position specifier for the text. Values of \code{1}, \code{2}, \code{3} and \code{4}, respectively indicate positions below, to the left of, above and to the right of the specified coordinates.}
#  \item{offset}{When `pos' is specified, this value gives the offset of the label from the specified coordinate in fractions of a character width.}
#  \item{xpd}{A logical value or @NA. If @FALSE, all plotting is clipped to the plot region, if @TRUE, all plotting is clipped to the figure region, and if @NA, all plotting is clipped to the device region.}
#  \item{...}{Common arguments accepted by most text labeling functions.
#   For more information see @see "graphics::text" and @see "graphics::par".}
# }
#
# \examples{
#   # Loads the file 'gpr123.gpr' located in the data directory:
#   gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#   raw <- getRawData(gpr);
#   ma <- getSignal(raw)
#
#   # Highlight and label the first 50 genes
#   idx <- 1:50
#   plot(ma)
#
#   # Highlight all spots in grid (1,2)
#   highlight(ma, idx, col="purple")
#
#   # Add bold faced labels using the labels specified by the layout
#   text(ma, idx, font=2)
# }
#
# @author
#
# \seealso{
#   @seemethod "plotXY".
#   @seemethod "plotSpatial".
#   @seemethod "highlight".
#   @see "graphics::text".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("text", "MicroarrayData", function(x, include=NULL, exclude=NULL, labels=NULL, col="black", offset=0.2, pos=4, xpd=TRUE, ...) {
  # To please R CMD check...
  this <- x;

  lastPlot <- Device$getPlotParameters();
  if (length(lastPlot) == 0)
    throw("To call points() of a MicroarrayData class one must first call one of the plot functions of a MicroarrayData class that creates a new plot must be called first.");
  fcn <- get(lastPlot$fcn, mode="function");

  if (is.null(labels)) {
    incl <- which(getInclude(this, include, exclude, slide=lastPlot$slide));
    layout <- getLayout(this);
    if (!is.null(layout))
      labels <- getID(layout, incl)
    else
      labels <- as.character(incl);
  }
  fcn(this, style="text", include=include, exclude=exclude, labels=labels, col=col, offset=offset, pos=pos, xpd=xpd, ...);
})


#########################################################################/**
# @RdocMethod lowessCurve
#
# @title "Draws one or more lowess curves through the data in last plot"
#
# \description{
#  @get "title".
#  Normally this function is only applicable to scatter plots.
# }
#
# @synopsis
#
# \arguments{
#  \item{include}{The indices of the spots that should be included. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are considered.
#   If @NULL all spots are considered.}
#  \item{exclude}{The indices of the spots that should be excluded. 
#   If it is instead a name of one or more flags, the spots which have been
#   flagged with these flags are excluded.
#   If @NULL no spots are excluded.}
#  \item{col}{The color(s) to be used for the plotted line(s). If @NULL
#   the global default line color will be used.}
#  \item{...}{Common arguments accepted by underlying plot functions.
#   For more information see @see "graphics::par" and @see "graphics::plot".}
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   plot(ma)
#   lowessCurve(ma, gridwise=TRUE)
# }
#
# @author
#
# \seealso{
#   @seemethod "text".
#   @see "graphics::points".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("lowessCurve", "MicroarrayData", function(this, include=NULL, exclude=NULL, col=NULL, f=0.3, ...) {
  lastPlot <- Device$getPlotParameters();
  if (length(lastPlot) == 0)
    throw("To call points() of a MicroarrayData class one must first call one of the plot functions of a MicroarrayData class that creates a new plot must be called first.");
  fcn <- get(lastPlot$fcn, mode="function");
  line <- fcn(this, style="lowessline", include=include, exclude=exclude, col=col, f=f, ...);
  invisible(line);
});


setMethodS3("lowessLine", "MicroarrayData", function(this, ...) {
  lowessCurve(this, ...)
}, private=TRUE, deprecated=TRUE)



############################################################################
# Color generators
############################################################################

setMethodS3("createColors", "MicroarrayData", function(this, x, y, type="MA", palette=NULL, A.range=c(0,16), M.range=c(-2,2)) {
  dim <- dim(x);
  
  if (!is.element(type, c("MA", "RG")))
    throw("Unknown type: ", type);
  
  if (is.null(palette)) palette <- "redgreen";
  if (!is.element(palette, c("redgreen", "blueyellow")))
    throw("Unknown palette: ", palette);

  if (is.null(A.range))
    A.range <- c(0,16);
  if (is.null(M.range))
    M.range <- c(-2,2);
  
  range <- matrix(c(A.range, M.range), nrow=2);

  if (type == "RG") {
    tmp <- 1/2*log(x*y, base=2);
    y <- log(x/y, base=2);
    x <- tmp;
    rm(tmp);
  }

  x <- matrix(c(x,y), ncol=2);
  x[is.na(x)] <- 0;
  
  if (palette == "redgreen") {
    dim.range <- matrix(c(0,1, Colors$GREEN.HUE,Colors$RED.HUE), nrow=2);
    col <- Colors$getHSV(x, x.range=range, dim=c("v", "h"), dim.range=dim.range);
  } else if (palette == "blueyellow") {
    dim.range <- matrix(c(0,1, Colors$YELLOW.HUE,Colors$BLUE.HUE), nrow=2);
    col <- Colors$getHSV(x, x.range=range, dim=c("v", "h"), dim.range=dim.range);
  }
  dim(col) <- dim;
  col;
}, static=TRUE, trial=TRUE)


# Generate grayscale colors by default.
setMethodS3("getColors", "MicroarrayData", function(this, what, slide=1, log=NULL, ...) {
#  if (!is.null(what) && !is.element(what, getFieldNames(this)))
#    throw("Argument 'what' is refering to an unknown field: ", what);
  what <- what[1];
  field <- this[[what]];

  view <- getView(this);
  if (is.null(view))
    view <- MicroarrayArray$DEFAULT.VIEW;
  if (view == MicroarrayArray$DEFAULT.VIEW)
    view <- getView(field);
  if (view == MicroarrayArray$SPOT.SLIDE.VIEW) {
    x <- getSpotSlideValues(field, slides=slide);
  } else if (view == MicroarrayArray$GENE.SLIDE.VIEW) {
    x <- getGeneSlideValues(field, slides=slide);
  } else if (view == MicroarrayArray$GENE.REPLICATE.SLIDE.VIEW) {
    x <- getGeneReplicateSlideValues(field, slides=slide);
    x <- unlist(x);
  } else {
    throw("Unsupported view.");
  }

  if (!is.null(log))
    x <- log(x, base=log);

  Colors$getGray(x);
})


#########################################################################/**
# @RdocMethod subplots
#
# @title "Prepare for a grid of subplots based on the number of slides"
#
# \description{
#  Creates a grid of subplots based on the number of slides this dataset
#  contains. 
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Any arguments that @see "R.graphics::subplots.Device" takes.}
# }
#
# \value{
#   Returns a @matrix with the shape of the grids and which contains the
#   subplot indices indicating the order that the subplots are plotted.
# }
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw)
#
#   subplots(ma)
#   for (k in seq(ma))
#     plot(ma, slide=k)
# }
#
# @author
#
# \seealso{
#   @see "graphics::layout".
#   @see "R.graphics::subplots.Device".
#   To set the distances between the subplots see graphical parameter 
#   \code{mar} in @see "graphics::par".
# }
#*/#########################################################################
setMethodS3("subplots", "MicroarrayData", function(this, ..., adjustMargins=TRUE, mar=c(2,2,1,1)) {
  res <- Device$subplots(nbrOfSlides(this), ...);
  if (adjustMargins) par(mar=mar);
  invisible(res);
})




############################################################################
# Special figure/plot annotations
############################################################################
if (is.null(getOption("com.braju.sma.plot.slide")))
  options(com.braju.sma.plot.slide=list(auto=TRUE, side=3, adj=1, cex=0.7, col="darkgray", line=0, outer=FALSE));

if (is.null(getOption("com.braju.sma.plot.gene")))
  options(com.braju.sma.plot.gene=list(auto=TRUE, side=3, adj=1, cex=0.7, col="darkgray", line=0, outer=FALSE, id="auto", name=FALSE));

setMethodS3("putSlide", "MicroarrayData", function(this, slide=NULL, auto=NULL, adj=NULL, cex=NULL, col=NULL, line=NULL, outer=FALSE, side=NULL, ...) {
  lastPlot <- Device$getPlotParameters();
  if (is.null(slide)) slide <- lastPlot$slide;
  if (is.null(slide)) {
    slideStr <- paste("All", nbrOfSlides(this));
  } else {
    slideStr <- paste("Slide:", paste(slide, collapse=",", sep=""));
    if (length(slide) == 1) {
      slideNames <- getSlideNames(this);
      if (!is.null(slideNames))
        slideStr <- paste(slideStr, " (", slideNames[slide], ")", sep="");
    }
  }

  o <- getOption("com.braju.sma.plot.slide");
  if (is.null(auto)) auto <- o$auto;
  if (is.null(auto) || !is.logical(auto) || !auto)
    return(invisible(FALSE));

  if (is.null(adj))   adj <- o$adj;
  if (is.null(adj))   adj <- 1;
  if (is.null(cex))   cex <- o$cex;
  if (is.null(cex))   cex <- 0.7;
  if (is.null(col))   col <- o$col;
  if (is.null(col))   col <- "darkgray";
  if (is.null(line))  line <- o$line;
  if (is.null(line))  line <- 0;
  if (is.null(outer)) outer <- o$outer;
  if (is.null(outer)) outer <- FALSE;
  if (is.null(side))  side <- o$side;
  if (is.null(side))  side <- 3;

  if (outer) line <- par("oma")[side] - line - 1;

  mtext(slideStr, side=side, adj=adj, cex=cex, col=col, line=line, outer=outer);
  return(invisible(TRUE));
}, trial=TRUE)


setMethodS3("putGene", "MicroarrayData", function(this, gene=NULL, auto=NULL, adj=NULL, cex=NULL, col=NULL, line=NULL, outer=FALSE, side=NULL, id=NULL, name=NULL, ...) {
  lastPlot <- Device$getPlotParameters();
  if (is.null(gene)) gene <- lastPlot$gene;
  if (is.null(gene)) gene <- paste("Unknown");
  geneStr <- "";

  o <- getOption("com.braju.sma.plot.gene");
  if (is.null(auto)) auto <- o$auto;
  if (is.null(auto) || !is.logical(auto) || !auto)
    return(invisible(FALSE));

  if (is.null(adj))   adj <- o$adj;
  if (is.null(adj))   adj <- 1;
  if (is.null(cex))   cex <- o$cex;
  if (is.null(cex))   cex <- 0.7;
  if (is.null(col))   col <- o$col;
  if (is.null(col))   col <- "darkgray";
  if (is.null(line))  line <- o$line;
  if (is.null(line))  line <- 0; 
  if (is.null(outer)) outer <- o$outer;
  if (is.null(outer)) outer <- FALSE;
  if (is.null(side))  side <- o$side;
  if (is.null(side))  side <- 3;
  if (is.null(id))    id <- o$id;
  if (is.null(id))    id <- FALSE;
  if (is.null(name))  name <- o$name;
  if (is.null(name))  name <- FALSE;

  if (id == TRUE && name == TRUE)
    warning("Only one of the arguments 'id' and 'name' can be TRUE. Discarding 'name'.");

  if (gene != "Unknown" && (id == TRUE || name == TRUE)) {
    layout <- getLayout(this);

    if (id && !hasIDs(layout)) id <- FALSE;
    if (name && !hasNames(layout)) name <- FALSE;

    if (hasReplicates(layout)) {
      reps <- getReplicates(layout);
      index <- getSpot(reps, gene)[1];
    } else {
      index <- gene;
    }
    if (id) {
      geneStr <- paste(" (", getID(layout, index), ")", sep="");
    } else if (name) {
      geneStr <- paste(" (", getName(layout, index), ")", sep="");
    }
  }

  geneStr <- paste("Gene: ", paste(gene, collapse=",", sep=""), geneStr, sep="");

  if (outer) line <- par("oma")[side] - line - 1;

  mtext(geneStr, side=side, adj=adj, cex=cex, col=col, line=line, outer=outer);
  return(invisible(TRUE));
}, trial=TRUE)




setMethodS3("plotGene", "MicroarrayData", function(this, what, gene=1, slides=NULL, pch=176, lty.within=1, lty.between=2, hline=median, ylog=NULL, new=TRUE, col="auto", ylab=what, xlab="gene replicate", ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  
  ys <- getGeneReplicateSlideValues(this[[what]], genes=gene, slides=slides);
  ys <- unlist(ys);
  nbrOfReplicates <- length(ys);

  if (!is.null(ylog) && ylog > 0)
    ys <- log(ys, base=ylog);

  if (!is.null(col) && col == "auto") {
    palette <- col;
    if (what %in% c("R", "Rb")) {
      if (max(ys, na.rm=TRUE) > 100)
        col <- Colors$getRed(log(ys,2), x.range=c(0,16))
      else
        col <- Colors$getRed(ys, x.range=c(0,16));
    } else if (what %in% c("G", "Gb")) {
      if (max(ys, na.rm=TRUE) > 100)
        col <- Colors$getGreen(log(ys,2), x.range=c(0,16))
      else
        col <- Colors$getGreen(ys, x.range=c(0,16));
    } else if (what == "M") {
      A <- rep(10, length(ys));
      col <- MicroarrayData$createColors(A, ys);
    } else if (what == "A") {
      col <- Colors$getGray(ys, x.range=c(0,16));
    }
  }

  Device$setPlotParameters(object=this, fcn="plotGene", what=what, gene=gene, ...);

  xs <- seq(ys);
  if (new == TRUE) {
    plot(xs, ys, pch=pch, col=col, xlab=xlab, ylab=ylab, ...)
    Device$putTimestamp();
    Device$putDataset();
    Device$putAuthor();
    putGene(this);
  } else {
    points(xs, ys, pch=pch, col=col, ...);
  }

  if (is.matrix(ys) && ncol(ys) > 1) {
    # Has within-slide replicates

    # Draw lines between within-slide replicates
		for (k in seq(ncol(ys))) {
 			xs <- seq(nrow(ys)*(k-1)+1, length=nbrOfReplicates);
 			lines(xs, ys[xs], lty=lty.within, ...);
    }

  	# Draw lines between between-slide replicates
 		for (k in seq(from=1, to=ncol(ys)-1)) {
 			xs <- seq(nrow(ys)*k, length=nbrOfReplicates);
 			lines(xs, ys[xs], lty=lty.between, ...);
  	}
  } else {
    # No within-slide replicates
  	# Draw lines between between-slide replicates
  	if (length(ys) > 1) {
 			xs <- seq(1, length=length(ys));
 			lines(xs, ys[xs], lty=lty.between, ...);
  	}
  }

  if (!is.null(hline)) 
    abline(h=hline(na.omit(ys)), lty=3, ...);

  if (!is.null(ys)) invisible(ys);
})


############################################################################
# HISTORY:
# 2005-05-04
# o Replace all image.matrix() with image270().
# 2004-10-22
# o BUG FIX: In R v2.0.0 lowessCurve(..., lwd=2) gives an error, because 
#   of duplicate matching of argument 'lwd'. Fixed by not drawing the gray
#   line with lwd=1.5*par("pwd") anymore, but with default width.
# 2004-05-22
# o Added Rdoc comments for plotDensity().
# 2004-02-17
# o Added plotDensity() dated 2003-11-10.
# 2004-01-14
# o Added the possibility to let 'xlim' and 'ylim' plotXY() be R expressions
#   too. This makes it possible to for instance calculate ylim <- 2*xlim
#   when xlim is known, which is sometimes calculated within plotXY(). This
#   now makes it possible to create M vs A plots where the M scale is twice
#   the A scale, which makes log R and log G ortogonal.
# 2003-07-30
# o Added Rdocs for plotSpatial3d().
# 2003-03-24
# o Renamed lowessLine() to lowessCurve() and made lowessLine() deprecated.
# 2002-12-21
# o BUG FIX: Plot annotation such as putting a timestamp or slide name was
#   also performed when highlight(), text() etc was called. Now it is only
#   done when the plot is created.
# 2002-12-11
# o Added more informative error messages if one tries to call for instance
#   points(<MAData>, ...) without calling plot(<MAData>, ...) first, but
#   only a standard plot(x,y) function. There are situations where the 
#   latter example should work, but for now the implementation does not
#   support it. We will work on this.
# 2002-12-06
# o Added support for specifying the groupBy argument using LayoutGroup
#   objects. The xlab is automatically constructed from the class name of
#   the groupBy object.
# 2002-12-05
# o Added argument 'groupBy=c("slide", "printtip")' for boxplot() to specify
#   how the "boxes" in the boxplot should represent. The boxes are now also
#   colored automatically if 'what' is "M" or "A".
# 2002-11-27
# o In boxplot() default value of 'gridwise' is now FALSE.
# 2002-11-07
# o Added more detailed source code comments about all the plot methods.
# o Added argument mode="function" to all get(lastPlot$fcn).
# o Added qqnorm().
# 2002-10-14
# o Added isFieldColorable().
# 2002-09-30
# o BUG in NextMethod()?: Can not use NextMethod() with plotPrintorder() or
#   plotDiporder(). Somehow not only what=what, but also slide=what. What
#   is going on???
# o Made plotPrintorder() public and added an reference to my manuscript.
# 2002-07-14
# * BUG FIX: Introduced a bug in plotSpatial() such that 'include' would give
#   an error if used.
# * Now points() for plotSpatial() also makes use of image().
# * Changed the default highlight() symbol for plotSpatial() to pch=20.
# 2002-06-30
# * plotSpatial() is now making use of image.matrix(), e.g. image(), instead
#   of plot(). Removed all anoying spacing between "spots" when the window
#   was resized.
# 2002-05-10
# * Added plotGene().
# 2002-05-07
# * Added plotDiporder().
# * Now createColors() also sets the dimensions of the returned color object
#   equal to the dimension of argument 'x'.
# 2002-05-01
# * Argument 'col' now sets the color of the lines in lowessCurve().
# * Now 'col' argument in plot() and plotSpatial() can be color numbers.
# * BUG FIX: Fixed the coloring of spots in plotPrintorder() when plotting
#   with include=from:to where from > 1.
# * Added points() and support for it in plotXY(), plotSpatial() and
#   plotPrintorder().
# 2002-04-21
# * Extracted MicroarrayData.PLOT.R
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

