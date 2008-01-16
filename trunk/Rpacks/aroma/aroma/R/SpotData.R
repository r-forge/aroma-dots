#########################################################################/**
# @RdocClass SpotData
#
# @title "The SpotData class"
#
# @synopsis
#
# \description{
#  @classhierarchy
#
#  Creates an SpotData object. If the data frame \code{data} is empty or
#  @NULL, the object will be empty.
# }
#
# \arguments{
#   \item{layout}{A @see "Layout" object specifying the spot layout of the
#    slides in this data set.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Details}{
#  A Spot file contains spot information for each spot on a single microarray slide. It consists of a header followed by a unspecified number of rows. The header contains 1+30 labels, and each row contains 31 fields. Each row corresponds to one spot. The fields are:\cr
#
#   \item{<NO NAME>}{row number}
#
#   \item{indexs}{spot number on slide. Range [0,N] in N.}
#
#   \item{grid.r}{grid row number. Range [1,GR] in Z+.}
#   \item{grid.c}{grid column number. Range [1,GC] in Z+.}
#   \item{spot.r}{spot row (within grid) number. Range [1,SR] in Z+.}
#   \item{spot.c}{spot column (within grid) number. Range [1,SC] in Z+.}
#
#   \item{area}{the number of foreground pixels. Range [0,MAXAREA] in N}
#
#   \item{Gmean}{the average of the foreground pixel values. Range [0,65535] in R}
#   \item{Gmedian}{the median of the foreground pixel values. Range [0,65535] in N}
#   \item{GIQR}{the inter quartile range (a robust measure of spread) of the logged foregroud pixel values. Range [0,16]+Inf in {R,NA}}
#   \item{Rmean}{the average of the foreground pixel values. Range [0,65535] in R}
#   \item{Rmedian}{the median of the foreground pixel values. Range [0,65535] in N}
#   \item{RIQR}{the inter quartile range (a robust measure of spread) of the logged foregroud pixel values. Range [0,16]+Inf in {R,NA}}
#
#   \item{bgGmean}{mean green background intesity. Range [0,65535] in R}
#   \item{bgGmed}{median green background intesity. Range [0,65535] in N}
#   \item{bgGSD}{standard deviation for the green background. Range [0,65535]+Inf in R}
#   \item{bgRmean}{mean red background intesity. Range [0,65535] in R}
#   \item{bgRmed}{median red background intesity. Range [0,65535] in N}
#   \item{bgRSD}{standard deviation for the red background. Range [0,65535]+Inf in R}
#
#   \item{valleyG}{the background intesity estimate from the local background valley method S.valley. Range [0,65535] in N}
#   \item{valleyR}{the background intesity estimate from the local background valley method S.valley. Range [0,65535] in N}
#
#   \item{morphG}{green background estimate using morphological opening (erosion-dilation). Range [0,65535] in N}
#   \item{morphG.erode}{green background estimate using morphological erosion. Range [0,65535] in N}
#   \item{morphG.close.open}{green background estimate using morphological closing-opening (dilation-erosion-dilation). Range [0,65535] in N}
#   \item{morphR}{red background estimate using morphological opening (erosion-dilation). Range [0,65535] in N}
#   \item{morphR.erode}{red background estimate using morphological erosion. Range [0,65535] in N}
#   \item{morphR.close.open}{red background estimate using morphological closing-opening (dilation-erosion-dilation). Range [0,65535] in N}
#
#   \item{logratio}{== log((Rmedian-morphR)/(Gmedian-morphG), base=2), i.e. \bold{Redundant}.}
#   \item{perimeter}{== 2*sqrt(pi*area/circularity), i.e. \bold{Redundant}.}
#   \item{circularity}{Shape of spot defined as 4*pi*area/perimeter**2.}
#   \item{badspot}{If the spot area is greater than product of the horizontal and the vertical average spot separations, equal to \code{1}, otherwise \code{0}.}
# }
#
# \section{About IQR}{
#   The interquartile range (IQR) is the distance between the 75\% 
#   @see "base::quantile" (percentile) and the 25\% quantile. 
#   In words, IQR is the range of the mid 50\%. Thus, no outliers are
#   included in the measure, which is why we say it is a robust measure.
#   For norammly distributed data \eqn{IQR = 1.35*\sigma}, where \eqn{\sigma}
#   is the standard deviation.
# }
# 
# \section{More on background estimates}{
#   The Spot software provides several different kinds of background estimates where three of them are 
#   based on morphological methods. For all of these methods, the signal selected to be the background
#   signal is the pixel value at the center of the spot \emph{after} applying the morphological transform
#   using a square mask with side 2.5 times the average distance between two spots.
#   The first and also the simpliest transform (\code{morph.erode}) performs a single \emph{erosion} step.
#   The second transform (\code{morph}) performs an \emph{opening}, which is an \emph{erosion} followed by
#   a \emph{dilution}. 
#   The third transform (\code{morph.close.open}) performs a \emph{closing} followed by an \emph{opening}, 
#   which is the same as doing a \emph{dilution}, then an \emph{erosion} and a \emph{dilution} again.
#   As the names of the steps indicate, an erosion makes the signal smaller and the dilution the signal
#   larger. Hence, background estimated based on these three methods can always be ordered as
#   \code{morph.erode <= morph <= morph.close.open}.
# }
#
# @author
#
# \references{
#   Spot Software package by CSIRO, Australia, 
#   \url{http://www.cmis.csiro.au/iap/spot.htm}
#
#   Spot: Description of Output, 2003
#   \url{http://www.cmis.csiro.au/iap/Spot/spotoutput.htm}
#
#   Y.H. Yang, M. Buckley, S. Dudoit, T. Speed, 
#   Comparison of methods for image analysis on cDNA microarray data, 
#   Tech. Report 584, Nov 2000.
#   \url{http://www.stat.berkeley.edu/users/terry/zarray/Html/image.html}
# }
#
# \examples{
#   spot <- SpotData$read("spot123.spot", path=system.file("data-ex", package="aroma"))
#
#   # Get the foreground and the background (and the layout)
#   raw <- getRawData(spot)
#
#   # The the background corrected data
#   ma <- getSignal(raw, bgSubtract=FALSE)
#
#   subplots(4, ncol=2)
#
#   # Plot R vs G with a lowess line through the data points
#   rg <- as.RGData(ma)
#   plot(rg)
#   lowessCurve(rg, lwd=2, gridwise=TRUE)
#
#   # Plot M vs A with a lowess line through the data points
#   plot(ma)
#   lowessCurve(ma, lwd=2, gridwise=TRUE)
#
#   # Plot spatial
#   plotSpatial(ma)
# }
#*/#########################################################################
setConstructorS3("SpotData", function(layout=NULL) {
  extend(MicroarrayData(layout=layout), "SpotData", 
    log.base = list()
  )
})


setMethodS3("append", "SpotData", function(this, other) {
  NextMethod("append");
  invisible(this);
})

setMethodS3("preferredLogBase", "SpotData", function(this, fields=getFieldNames(this)) {
  logbase <- list(indexs=0, grid.r=0, grid.c=0, spot.r=0, spot.c=0, area=0, Gmean=2, Gmedian=2, GIQR=2, Rmean=2, Rmedian=2, RIQR=2, bgGmean=2, bgGmed=2, bgGSD=2, bgRmean=2, bgRmed=2, bgRSD=2, valleyG=2, valleyR=2, morphG=2, morphG.erode=2, morphG.close.open=2, morphR=2, morphR.erode=2, morphR.close.open=2, logratio=0, perimeter=0, circularity=0, badspot=0);
  unlist(logbase[fields]);
}, trial=TRUE)


setMethodS3("logPreferred", "SpotData", function(this) {
  log.base <- preferredLogBase(this);
  field <- getFieldNames(this);
  for (k in seq(length(log.base))) {
    name <- field[k];
    if (log.base[k] > 0) {
      this[[name]] <- log(this[[name]], base=log.base[k]);
      this$log.base[[name]] <- log.base[k];
    }
  }
  invisible(this);
}, trial=TRUE)


# Is this needed? /HB 2005-02-28
#if (!exists("log.default"))
#  log.default <- log;

setMethodS3("log", "SpotData", function(x, field, base=2, ...) {
  # To please R CMD check...
  this <- x;

  base <- rep(base, length.out=length(field));
  for (k in 1:length(field)) {
    name <- field[k];
    if (base[k] > 0) {
      this[[name]] <- log(this[[name]], base=base[k]);
      this$log.base[[name]] <- base[k];
    }
  }
  invisible(this);
})




#########################################################################/**
# @RdocMethod read
#
# @title "Reads Spot Results Data file(s)"
#
# \description{
#  Reads file(s) that are of the Spot file format (Spot Results Format).
# }
#
# @synopsis
#
# \arguments{
#   \item{filenames}{The filename(s) of the Spot file(s) to be read.}
#   \item{path}{The path(s) to the Spot file(s).}
#   \item{sep}{Character that separates the columns of data.
#     If \code{"auto"} a best guess is made from information in the header.}
#   \item{skip}{Number of lines to skip \emph{before} reading the header.
#     If \code{"auto"}, the number of lines to skip before the header will
#     be searched for automatically.}
#   \item{verbose}{If @TRUE, information will printed out during
#                  the file reading/parsing.}
# }
#
# \value{
#  Returns a @see "SpotData" object containing one or several slides.
# }
#
# @author
#
# \examples{\dontrun{# For an example see help(SpotData).}}
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("readOneFile", "SpotData", function(this, filename, path=NULL, sep="auto", skip="auto", verbose=TRUE) {
  filename <- Arguments$getReadablePathname(filename, path);  

  if (verbose) cat("Reading file ", filename, "...", sep="");

  if (skip == "auto" || sep == "auto") {
    # Support gzip'ed files too.
    lines <- readLines(gzfile(filename), n=20);
    if (skip == "auto") {
      skip <- grep("indexs", lines);
      if (length(skip) == 0)
        skip <- 0
      else
        skip <- skip - 1;
    }
    if (sep == "auto") {
      header <- lines[skip+1];
      len0 <- nchar(header);
      lenCommas <- len0 - nchar(gsub(",","", header));
      lenTABs   <- len0 - nchar(gsub("\t","", header));
      if (lenCommas > lenTABs)
         sep <- ","
      else
         sep <- "\t";
    }
  }

  # Support gzip'ed files too.
  if (regexpr("[.]gz$", filename) != -1) {
    tmpname <- tempfile();
    n <- gunzip(filename, tmpname);
    filename <- tmpname;
    on.exit(file.remove(tmpname));
  }

  # Read the data with the header
  # The header/data *might* be quoted with " or ', but could be unquoted too.
  # Support gzip'ed files too.
  df <- read.table(file=filename, header=TRUE, sep=sep, quote="\"'", skip=skip, strip.white=TRUE, blank.lines.skip=TRUE);

  # Rename the field names to standard field names.
  colnames(df) <- SpotData$renameFields(colnames(df));
  
  # Get the layout.
  gridfields <- c("grid.r", "grid.c", "spot.r", "spot.c");
  layout <- SpotData$extractLayout(df[gridfields]);

  # Create 
  spot <- SpotData(layout=layout)
  spot$.fieldNames <- names(df);
  for (field in names(df)) {
    spot[[field]] <- as.matrix(df[[field]]);
    df[[field]] <- NULL; # Save memory
  }

  if (verbose) cat("ok\n", sep="");
  
  rm(df); gc(); # To minimize memory usage!
  
  spot;
}, protected=TRUE, static=TRUE);



#########################################################################/**
# @RdocMethod write
#
# @title "Write a SpotData object to file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file to be written.}
#   \item{path}{The path to the file.}
#   \item{slide}{Index of slide to be saved.}
#   \item{...}{Other arguments accepted by subclasses or which are passed
#     to \code{write.table}.}
# }
#
# \value{Returns nothing.}
#
# @author
#
# \examples{
#   spot <- SpotData$read("spot123.spot", path=system.file("data-ex", package="aroma"))
#
#   # Write the SpotData object to a temporary file.
#   filename <- paste(tempfile("SpotData"), ".dat", sep="")
#   write(spot, filename)
#
#   spot2 <- SpotData$read(filename)
#   print(equals(spot, spot2))  # TRUE
#
#   unlink(filename)
# }
#
# \seealso{
#   To read one or more SpotData files at once
#   see @see "SpotData.read".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("write", "SpotData", function(this, filename, path=NULL, slide=NULL, overwrite=FALSE, row.names=FALSE, ..., verbose=FALSE) {
  filename <- Arguments$getWritablePathname(filename, path, mustNotExist=!overwrite);  

  slide <- validateArgumentSlide(this, slide=slide);

  df <- extract(this, slides=slide);
  write.table(df, file=filename, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t", append=FALSE, ...);
});




#########################################################################/**
# @RdocMethod read
#
# @title "Reads several Spot files into a SpotData object"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{A @vector of filenames. Either \code{pattern} or \code{filename} must be specified.}
#   \item{path}{A string (or an optional @vector of paths if \code{filename} is specified) to the files.}
#   \item{pattern}{A pattern string for matching filenames. Either \code{pattern} or \code{filename} must be specified.}
#   \item{verbose}{If @TRUE, information will printed out during
#                  the reading of the file.}
# }
#
# \value{Returns a @see "SpotData" object.}
#
# @author
#
# \examples{\dontrun{# For an example see help(SpotData).}}
#
# \seealso{
#   For pattern formats see @see "base::list.files".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("read", "SpotData", function(this, filename=NULL, path=NULL, pattern=NULL, verbose=TRUE, ...) {
  if (is.null(filename) && is.null(pattern))
    throw("Either 'filename' or 'pattern' must be specified.");
  if (!is.null(filename) && !is.null(pattern))
    throw("Only one of 'filename' and 'pattern' can be specified.");

  if (is.list(path)) {
    path <- lapply(path, FUN=as.character);
  } else {
    path <- sapply(path, FUN=as.character);
  }

  # Make sure 'path' is a vector; note that sapply(NULL, ...) gives list()!
  path <- unlist(path);

  if (!is.null(pattern)) {
    # Remove '/' at the end since Windows system does not support this as a path.
    # I have discussed this with [R] developers, but they think it should be like
    # that, but it is not cross plattform safe.
    if ((pos <- regexpr(".*/$", path)) != -1)
      path <- substring(path, 1, pos-1);
    if (is.null(path) || path == "") path0 = "." else path0 <- path;
    filename <- list.files(path=path0, pattern=pattern);
    if (length(filename) == 0)
      return(list());
    if (verbose) cat("Loading ", length(filename), " files:\n", sep="");
  }
  
  res <- list();

  # Support both path as a single string or as a vektor of strings.
  path <- rep(path, length.out=length(filename));

  res <- NULL;
  for (k in 1:length(filename)) {
    gc(); # Call the garbage collector in case we are running low in memory.
    tmp <- SpotData$readOneFile(filename=filename[k], path=path[k], verbose=verbose, ...);
    slidename <- basename(filename[k]);
    setSlideNames(tmp, slidename);
    if (is.null(res))
      res <- tmp
    else {
      append(res, tmp);
    }
    rm(tmp);
  }

  gc();

  res;
}, static=TRUE);


# Kept for backward compatibility.
setMethodS3("readAll", "SpotData", function(...) {
  SpotData$read(...);
}, private=TRUE, static=TRUE, deprecated=TRUE)




setMethodS3("extractLayout", "SpotData", function(this, griddata) {
  ngrid.r <- max(griddata$grid.r);
  ngrid.c <- max(griddata$grid.c);
  nspot.r <- max(griddata$spot.r);
  nspot.c <- max(griddata$spot.c);
  Layout(ngrid.r, ngrid.c, nspot.r, nspot.c)
}, static=TRUE);



setMethodS3("renameFields", "SpotData", function(this, header) {
  map <- c("grid_r"="grid.r", "grid_c"="grid.c", "spot_r"="spot.r", "spot_c"="spot.c");
  map <- c(map, "shape"="circularity");
  map <- c(map, "cy3mean"="Gmean", "cy3IQR"="GIQR", "cy5mean"="Rmean", "cy5IQR"="RIQR", "bg3mean"="bgGmean", "bg3med"="bgGmed", "bg3SD"="bgGSD", "bg5mean"="bgRmean", "bg5med"="bgRmed", "bg5SD"="bgRSD", "valley3"="valleyG", "valley5"="valleyR");
  
  # Uncertain about these guys /HB 020910
  map <- c(map, "bgGIQR"="bgGSD", "bgRIQR"="bgRSD");
  map <- c(map, "cy3bg2"="morphG.erode", "cy5bg2"="morphR.erode", "cy3bg3"="morphG", "cy5bg3"="morphR");

  # There was one very special version that has to be handles like this:
  if (is.element("slide.c", header)) {
    map <- c(map, "slide.r"="grid.r", "slide.c"="grid.c", "grid.r"="spot.r", "grid.c"="spot.c");
  }

  for (k in seq(along=header)) {
    name <- header[k];
    newName <- map[name];
    if (!is.na(newName) && !is.null(newName))
      header[k] <- newName;
  }
  
  return(header);
}, private=TRUE, static=TRUE);





#########################################################################/**
# @RdocMethod getRawData
#
# @title "Gets the raw intensites from the SpotData structure"
#
# \description{
#  Extracts the red and green spot intensitites (both foreground and background)
#  from the SpotData object and returns a @see "RawData" object.
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{Specifying which slides to be extracted. If @NULL, 
#     all slides are considered.}
#   \item{fg}{If \code{"mean"}, the mean foreground intensities are returned.
#     If \code{"median"}, the median foreground intensities are returned.}
#   \item{bg}{If \code{"mean"}, the mean background intensities are returned. 
#     If \code{"median"}, the median background intensities are returned. To 
#     get the mean of the median of the background intensites in the four 
#     diamond shaped ares around the spots, use the value \code{"valley"}. 
#     There are also three morphological background estimates; \code{morph},
#     \code{morph.erode}, and \code{morph.close.open}.}
# }
#
# \value{
#   Returns a @see "RawData" object containing the specified slides.
# }
#
# \details{
#   The R and Rb channels will come from the *R* fields, and
#   the G and Gb channels will come from the *G* fields.
#   To swap the channels use dyeSwap().
# }
#
# \examples{\dontrun{# For an example see help(SpotData).}}
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getRawData", "SpotData", function(this, slides=NULL, fgs="auto", bgs="auto") {
  slides <- validateArgumentSlides(this, slides=slides);
  fieldNames <- getFieldNames(this);
  
  if (fgs == "auto") {
    if (is.element("Rmedian", fieldNames)) {
      fgs <- c("Rmedian", "Gmedian");
    } else if (is.element("Rmean", fieldNames)) {
      fgs <- c("Rmean", "Gmean");
    }
  }
  if (bgs == "auto") {
    if (is.element("morphR.close.open", fieldNames)) {
      bgs <- c("morphR.close.open", "morphG.close.open");
    } else if (is.element("morphR", fieldNames)) {
      bgs <- c("morphR", "morphG");
    } else {
      bgs <- c("morphR.erode", "morphG.erode");
    }
  }
  
  nas <- which(is.na(match(c(fgs,bgs), fieldNames)));
  if (length(nas) > 0) {
    throw("Strange, some fields do not exists in the SpotData object. Please, report this error to the author of the package: ", c(fgs,bgs)[nas]);
  }

  R  <- this[[fgs[1]]][,slides];
  G  <- this[[fgs[2]]][,slides];
  Rb <- this[[bgs[1]]][,slides];
  Gb <- this[[bgs[2]]][,slides];

  RawData(R=R, G=G, Rb=Rb, Gb=Gb, layout=getLayout(this), 
                                        extras=this$.extras)
})


setMethodS3("as.RawData", "SpotData", function(this, ...) {
  getRawData(this, ...);
})




#########################################################################/**
# @RdocMethod plotSpatial
#
# @title "Creates a spatial plot of a slide"
#
# \description{
#   @get "title".
#  Note that older versions of the
#  Spot software did not generate result file containing information about
#  the spatial coordinates of the spots. If this is the case, a standard
#  spatial plot is generated.
# }
#
# @synopsis
#
# \arguments{
#   \item{slide}{The slide to be plotted.}
#   \item{pch}{The spot symbol. Default is \code{20} (solid disk).}
#   \item{yaxt, ...}{Parameters accepted by \code{plot}.}
# }
#
# \value{Returns nothing.}
#
# @author
#
# \examples{\dontrun{For an example see help(SpotData).}}
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("plotSpatial", "SpotData", function(this, what=NULL, slide=1, include=NULL, exclude=NULL, pch=20, yaxt=NULL, cex=NULL, col="auto", palette="redgreen", A.range=c(0,16), M.range=c(-1,1), look=c("real", "classic"), style=NULL, ...) {
  look = match.arg(look);

  cmd <- NULL;
  if (!is.null(style) && is.element(style, c("points", "highlight", "text"))) {
    cmd <- style;
    lastPlot <- Device$getPlotParameters();
    what <- lastPlot$what;
    slide <- lastPlot$slide;
    look <- lastPlot$look;
    if (is.null(look))
      look <- "classic";
  }

  if (look == "classic") {
    if (is.null(what))
      what <- "logratio";
    return(plotSpatial.MicroarrayData(this, what=what, slide=slide, include=include, exclude=exclude, col=col, style=style, ...));
  }

  xField <- "nominal.center.col";
  yField <- "nominal.center.row";
  if (!hasField(this, xField) || !hasField(this, yField)) {
    warning(paste("Could not find information about spot coordinates. Using 'classic' look instead: ", xField, ", ", yField, sep=""));
    if (is.null(what))
      what <- "logratio";
    return(plotSpatial.MicroarrayData(this, what=what, slide=slide, include=include, exclude=exclude, col=col, style=style, ...));
  }
  
  if (length(slide) != 1 || !is.numeric(slide))
    throw("Argument 'slide' must be an integer: ", slide, collapse=", ");
  
  if (!is.null(what) && !is.element(what, getFieldNames(this)))
    throw("Argument 'what' is not refering to a known field: ", what);

  if (is.null(col) || col == "auto") {
    if (is.null(what)) {
      R <- this[["Rmedian"]][,slide];
      G <- this[["Gmedian"]][,slide];
      col <- MicroarrayData$createColors(R,G, type="RG", palette=palette, A.range=A.range, M.range=M.range);
    } else {
      col <- getColors(this, what, slide=slide);
    }
  }

  setView(this, MicroarrayArray$SPOT.SLIDE.VIEW);
  include <- which(getInclude(this, include, exclude, slide = slide));
  colWhat <- what;
  if (length(col) == 1) {
    color <- col;
    col <- rep(NA, length.out=nbrOfSpots(this));
    if (!is.numeric(color) && substring(color, 1, 1) != "#" &&
        !is.element(color, colors())) {
        col[include] <- getColors(this, what=colWhat, slide=slide,
                                  include=include, palette=color, log=log);
    }
    else {
        col[include] <- color;
    }
  }  

  xy <- getSpotPosition(this);
  # Upper *left* corner is (0,0)
  x <- xy$x; y <- -xy$y; # Will also be used as xlab and ylab.
  
  if (missing(yaxt))
    yaxt0 <- "n";

  if (is.null(cmd)) {
    plot(x,y, pch=pch, yaxt=yaxt0, cex=cex, col=col, ...);
  
    if (missing(yaxt)) {
      # Add the y-axis with the correct ticks
      yaxp <- par("yaxp");
      yaxis <- seq(yaxp[1],yaxp[2], length=yaxp[3]+1);
      axis(2, at=yaxis, labels=-yaxis);
    }
  } else if (cmd == "highlight") {
    if (length(pch) == 1 && pch == "auto")
      pch <- 20;
    col <- col[include];
    x <- x[include];
    y <- y[include];
    points(x,y, cex=cex, col=col, pch=pch, ...)
  }

  Device$setPlotParameters(object=this, fcn="plotSpatial", what=what, slide=slide, look="real");
}) # plotSpatial()


setMethodS3("plotSpatial3d", "SpotData", function(this, field="logratio", ...) {
  NextMethod("plotSpatial3d", this, field=field, ...);
})





setMethodS3("normalizeGenewise", "SpotData", function(this, fields=NULL, bias=0, scale=1, ...) {
  if (is.null(fields)) {
    fields <- getFieldNames(this);
    exclude <- c("indexs", "grid.r", "grid.c", "spot.r", "spot.c", "badspot");
    fields <- setdiff(fields, exclude);
  }
  NextMethod("normalizeGenewise", this, fields=fields, bias=bias, scale=scale, ...);
})


setMethodS3("getBackground", "SpotData", function(this, which=c("morph.close.open", "morph", "morph.erode", "valley", "bgmed", "bgmean")) {
  which <- match.arg(which);

  if (which == "morph.erode") {
    fields <- c("morphR.erode", "morphG.erode");
  } else if (which == "morph") {
    fields <- c("morphR", "morphG");
  } else if (which == "morph.close.open") {
    fields <- c("morphR.close.open", "morphG.close.open");
  } else if (which == "valley") {
    fields <- c("valleyR", "valleyG");
  }

  # Assert that the fields do really exist.
  if (!all(fields %in% getFields(this))) {
    throw("The background estimates ", paste(fields, collapse=" and "), 
          " is not part of this ", data.class(this), " object.");
  }

  RGData(R=this[[fields[1]]], G=this[[fields[2]]], layout=getLayout(this));
})


setMethodS3("getForeground", "SpotData", function(this, which=c("median","mean")) {
  which <- match.arg(which);

  if (which == "mean") {
    fields <- c("Rmean", "Gmean");
  } else if (which == "median") {
    fields <- c("Rmedian", "Gmedian");
  }

  # Assert that the fields do really exist.
  if (!all(fields %in% getFields(this))) {
    throw("The background estimates ", paste(fields, collapse=" and "), 
          " is not part of this ", data.class(this), " object.");
  }

  RGData(R=this[[fields[1]]], G=this[[fields[2]]], layout=getLayout(this));
})



setMethodS3("getArea", "SpotData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  as.matrix(this[["area"]][include,slides]);
})


setMethodS3("getCircularity", "SpotData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  as.matrix(this[["circularity"]][include,slides]);
})


setMethodS3("getSNR", "SpotData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  snr.R <- this[["s2n5"]][include,slides];
  snr.G <- this[["s2n3"]][include,slides];
  snr <- sqrt(snr.R * snr.G);
  as.matrix(snr);
})


setMethodS3("getPerimeter", "SpotData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  as.matrix(this[["perimeter"]][include,slides]);
})


setMethodS3("getBgArea", "SpotData", function(this, slides=NULL, include=NULL, ...) {
  throw("Sorry, but data from Spot does not contain information about the number of pixels used in the background estimation.");

  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));
})


#########################################################################/**
# @RdocMethod getForegroundSD
#
# @title "Gets (an approximation of) the standard deviation of the foreground pixels"
#
# @synopsis
#
# \description{
#  @get "title". 
#  Since the Spot software does not return estimates of the standard
#  deviation of the foreground pixels, but the interquartile range (IQR) of
#  their logarithm instead, the value returned by this method is only an
#  approximation of the standard deviation. 
#  \emph{It is based on the assumption that the noise is symmetric on the
#  non-logaritmic scale}.
# }
#
# \value{
#   Returns a list of matrices that contain approximations of the estimated
#   standard deviation of the pixels in the foreground region of the spots.
# }
#
# \details{
#  By assuming that the noise is symmetric on the non-logaritmic scale we
#  can first derive the IQR on the non-logaritmic scale since we know the
#  median pixel value too. Let the non-logaritmic pixel value be \eqn{X} and
#  the let \eqn{x=log[2](X)}. More over, let \eqn{X[0.50]} be the median
#  value of all \eqn{X}'s, \eqn{X[0.25]} and \eqn{X[0.75]} be the 0.25 and
#  0.75 quantile, respectively. \eqn{x[0.25], x[0.50], x[0.75]} are defined
#  analogously. Assuming symmetric noise on the non-log scale means that
#  \eqn{X[0.50]-X[0.25] == X[0.75]-X[0.50]} is true.
#  Let \eqn{dX = X[0.75]-X[0.50]}. Given the IQR or the log pixel values,
#  i.e. \eqn{xIQR=x[0.75]-x[0.25]}, one can derive 
#  \deqn{
#    dX = \frac{2^{x_{IQR}}-1}{2^{x_{IQR}}+1} \cdot X_{0.50}
#  }{
#    dX = (2^xIQR-1/(2^xIQR+1)*X[0.50]
#  } 
#  Hence, the IQR of the non-logarithmic pixel values can then be
#  approximated by \eqn{2*dX}. If the noise is not only symmetric, 
#  but also Gaussian, the standard deviation is \eqn{1.35*2*dX}, 
#  which is the value returned.
# }
#
# \examples{
#   spot <- SpotData$read("spot123.spot", path=system.file("data-ex", package="aroma"))
#   raw <- getRawData(spot)
#   sd <- getForegroundSD(spot)
#   raw$RSD <- sd$RSD; raw$GSD <- sd$GSD;
#   subplots(4)
#   plot(raw, "RSDvsR", col="red")
#   plot(raw, "GSDvsG", col="green")  
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################  
setMethodS3("getForegroundSD", "SpotData", function(this) {
#  warning("SpotData does not contain estimates of the standard deviation of the foreground pixels, but the interquartile range (IQR) of the logarithm of them. The value returned by getForegroundSD() is an estimate based on the assumption that the noise is symmetric on the non-logaritmic scale.");

  res <- list();
  for (ch in c("R", "G")) {
    lxIQR <- this[[paste(ch, "IQR", sep="")]];
    xIQR <- 2^lxIQR;
    rm(lxIQR);
    factor <- (xIQR-1)/(xIQR+1);
    rm(xIQR);
    xMedian <- paste(ch, "median", sep="");
    if (!is.element(xMedian, names(this))) {
      msg <- paste("The field ", xMedian, " was not found. ", sep="");
      xMedian <- paste(ch, "mean", sep="");
      msg <- paste("Approximated by ", xMedian, " instead.", sep="");
      warning(msg);
    }
    res[[paste(ch, "SD", sep="")]] <- 1.35 * 2 * factor * this[[xMedian]];
  }
  res;
}) # getForegroundSD()


#########################################################################/**
# @RdocMethod getForegroundSE
#
# @title "Gets the standard error of the foreground pixels"
#
# @synopsis
#
# \description{
#  @get "title". 
# }
#
# \value{
#   Returns a list of matrices that contain the standard error
#   of the pixels in the foreground region of the spots.
# }
#
# \details{
#   The standard error returns the standard deviation divided by the area.
# }
#
# \examples{
#   gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#   raw <- getRawData(gpr)
#   sd <- getForegroundSE(gpr)
#   raw$RSE <- sd$RSE; raw$GSE <- sd$GSE;
#   subplots(4)
#   plot(raw, "RSEvsR", col="red")
#   plot(raw, "GSEvsG", col="green")
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################  
setMethodS3("getForegroundSE", "SpotData", function(this) {
  sd <- getForegroundSD(this);
  n <- getArea(this);
  sd <- lapply(sd, FUN=function(x) x / n);
  names(sd) <- gsub("SD", "SE", names(sd));
  sd;
})



#########################################################################/**
# @RdocMethod getSpotPosition
#
# @title "Gets physical positions of the spots"
#
# \description{
#  Gets physical positions (in pixels) of the spots on one or several
#  slides. 
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{Specifying which for slides the spot positions should
#     be extracted. If @NULL, all slides are considered.}
#   \item{index}{The spots for which the position is returned.
#     If @NULL all spots are considered.}
# }
#
# \value{Returns a @see "SpotPosition" object containing the
#   positions of the spots on the specified slides.}
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getSpotPosition", "SpotData", function(this, slides=NULL, index=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (any(index < 1) || any(index > nbrOfSpots(this)))
    throw("Argument 'index' is out of range.");
  if (is.null(index))
    index <- 1:nbrOfSpots(this);

  xField <- "nominal.center.col";
  yField <- "nominal.center.row";
  if (!hasField(this, xField) || !hasField(this, yField)) {
    throw("This ", data.class(this), " object is missing the fields '", xField, "' and/or '", yField, "', which specifies the physical position of the spots.");
  }

  x <- this[[xField]][index,slides];
  y <- this[[yField]][index,slides];
  SpotPosition(x=x, y=y);
})


############################################################################
#
#   A r r a y   d i m e n s i o n s
#
############################################################################
setMethodS3("getArrayLeft", "SpotData", function(this, slide=NULL) {
  slide <- validateArgumentSlide(this, slide=slide);
  getLeftEdge(getSpotPosition(this, slide=slide));
})

setMethodS3("getArrayRight", "SpotData", function(this, slide=NULL) {
  slide <- validateArgumentSlide(this, slide=slide);
  getRightEdge(getSpotPosition(this, slide=slide));
})

setMethodS3("getArrayTop", "SpotData", function(this, slide=NULL) {
  slide <- validateArgumentSlide(this, slide=slide);
  getTopEdge(getSpotPosition(this, slide=slide));
})

setMethodS3("getArrayBottom", "SpotData", function(this, slide=NULL) {
  slide <- validateArgumentSlide(this, slide=slide);
  getBottomEdge(getSpotPosition(this, slide=slide));
})

setMethodS3("getArrayWidth", "SpotData", function(this, slide=NULL) {
  getMaxWidth(getSpotPosition(this, slide=slide));
})

setMethodS3("getArrayHeight", "SpotData", function(this, slide=NULL) {
  getMaxHeight(getSpotPosition(this, slide=slide));
})

setMethodS3("getArrayAspectRatio", "SpotData", function(this, slide=NULL) {
  getAspectRatio(getSpotPosition(this, slide=slide));
})


############################################################################
# HISTORY:
# 2008-01-15
# o BUG FIX: Added a 'colWhat <- col' to plotSpatial() of SpotData.
# 2005-10-21
# o Replace 'overwrite' arguments with 'mustNotExist' in calls to Arguments. 
# 2005-07-20
# o BUG FIX: The internal 'path' of read() was a list, not a vector.
# 2005-07-19
# o Replaced all path="" arguments to path=NULL.
# 2005-06-11
# o Making use of Arguments in R.utils.
# 2005-05-03
# o Updated regular expressions.
# 2005-02-28
# o TYPO: All getArrayNNN() methods were defined for the GenePixData class
#   and not the SpotData class.
# 2004-??-??
# o Added support for File object as paths to read().
# 2004-08-15
# o Renamed as.RawData() to getRawData().
# 2004-03-09
# o BUG FIX: write() was mistakenly accepting more than once slide at the
#   time. It was a typo and now an error is thrown is more slides are given.
# 2003-10-30
# o Added methods for retrieving the array dimensions.
# o Added getSpotPosition().
# 2003-10-26
# o When trying to read a *.gz it is extracted to a temporary file, which is
#   read instead. Thus, *.spot.gz files are supported now.
# 2003-10-20
# o Added gzfile() to automatically identify the size of the header. Still
#   have to figure out how to make read.table() to support gzfile().
# 2003-10-19
# o getBackground() now returns morphological close open by default. Before
#   morph.erode() was returned. getForeground() returns *median* estimate
#   instead of mean. as.RawData() already had these.
# 2003-06-15
# o Added getForegroundSE() and getForegroundSD().
# o Updated the Rdoc about 'badspot' and added information about IQR for
#   non-statisticians.
# 2003-05-13
# o Update the help about 'logratio'; it is not calculated based on the 
#   'Rmean' but the 'Rmedian'. It is still redundant though.
# 2003-04-12
# o Added getBackground() and getForeground() for convienience.
# 2003-02-25
# o BUG FIX: Updated usage of read().
# 2003-01-31
# o Updated the Rdoc comments for SpotData to explain the three different
#   morphological background estimates much better.
# 2002-12-24
# o Added plotSpatial3d().
# o Updated the plotSpatial() method to generate a classical spatial plot
#   in case there is no information about the spot coordinates.
# 2002-12-21
# o When reading slides from files each slide is now named as the filename.
# 2002-12-11
# o The background estimates used by as.RawData() are now (in order)
#   "morph*", "morph*.close.open" and "morph*.erode". Before it was only
#   the last two that were used.
# o BUG FIX: read() would not extract the correct headers of a Spot files
#   if they where quouted. 
# 2002-12-05
# o BUG FIX: The get*() methods did not return a matrix if there was only
#   one slide.
# 2002-10-11
# o BUG FIX: morph.close.open was not always used by as.RawData().
# 2002-09-24
# o Changed the attribute 'path="."' to 'path=""' in read().
# o Updated the Rdoc comments for as.RawData().
# 2002-09-21
# o read() can now read one or several files specified by names or by a
#   pattern. This is identical to readAll(), which is now just calling
#   read() for backward compatibilities.
# 2002-09-12
# o BUG FIX: Some spot files do include a first row numbering column and
#   others do not. Since I before tried to use scan() to read the data I had
#   to make different tests for the above. However, read.table() deals with
#   such situations automatically => removed readInternal().
# 2002-09-10 [v0.49]
# o Added support for automatically detecting 'sep' and 'skip' in read().
# o Removed the private 'version' field. Not needed anymore.
# o Removed getVersion(), identifyVersion() and standardizeHeaders().
#   Now the read() methods is much more forgiving and accept new unknown
#   fields too.
# o Added renameFields().
# 2002-08-20
# o Replace 'append(super(this),other)' with NextMethod("append") in
#   method append().
# 2002-05-03
# o Added getArea(), getCircularity() and getSNR(). Also added Spot specific
#   getPerimeter(). The "cool" thing is that these can now be accessed as
#   spot$area, spot$circularity, spot$SNR etc. Assignment is still not
#   supported, i.e. READ ONLY.
# 2002-04-21
# * Added trial version of normalizeGenewise().
# 2002-04-05
# * Removed getSpotPositions(). There is a getSpotPosition() in the super
#   class which only make use of the Layout object, but it is better than
#   nothing.
# 2002-04-03
# * read() and write() now supports File objects.
# 2002-03-25
# * Added "dummies" for getSpotPositions() and plotSpatial().
# 2002-02-28
# * Made class implements interface Morpable.
# 2002-02-26
# * Added the write() method.
# * Removed as.character() from this class since there is now a generic
#   as.character() in the class MicroarrayData.
# * Modified code to make use of setMethodS3().
# 2001-11-17
# * Updated readAll() to also include pattern matching.
# 2001-11-12
# * Added readAll.SpotData().
# 2001-09-20
# * Added support for two old comma separated formats of Spot data.
# 2001-08-09
# * Modified readInternal to basically just read the file into a data frame.
#   The constructor the takes care of the headers, removing redundant 
#   fields are checking for version format. The constructor now accepts the
#   a data frame such as mouse1 in data( MouseArray).
# 2001-07-12
# * Bug fix: Forgot to do quote="" in all scan() and read.table() calls.
#   Trying to read files with cells including "'" would give an error.
# 2001-07-11
# * Renamed the class from SpotResultsData to SpotData.
# 2001-06-30
# * Renamed getRawData() to as.RawData(). Update to also include the layout.
# * Update some of the Rdoc comments.
# 2001-05-14
# * The .Internal type.convert() produced factors when data contained "Inf",
#   which resulted in way to big data.frames! A work around was necessary. 
#   See solution in readAll().
# * Added getInternalReferences() for improving gco() performance.
# 2001-05-13
# * Added append() and nbrOfSlides().
# * Added getLayout(), nbrOfSpots(), nbrOfFields(). Made all fields private.
# 2001-04-13
# * Created from GenePixResultsData.
############################################################################
