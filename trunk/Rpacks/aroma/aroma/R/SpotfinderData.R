#########################################################################/**
# @RdocClass SpotfinderData
#
# @title "The SpotfinderData class"
#
# @synopsis
#
# \description{
#  @classhierarchy
#
#  Creates an SpotfinderData object. If the data frame \code{data} is empty or
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
#  A Spotfinder/TAV (TIGR Array Viewer) file is a headerless tab-delimited 
#  file with file extension \code{*.tav}, which the following contains 8 or 
#  17 columns:\cr
#
#  The first eigth column do always exist:\cr
#   \item{1}{spot row number. Range [0,ROWS] in Z+.}
#   \item{2}{spot column number. Range [0,COLS] in Z+.}
#   \item{3}{spot metarow/subrow number. Range [0,MROWS] in Z+.}
#   \item{4}{spot metacolumn/subcolumn number. Range [0,MCOLS] in Z+.}
#   \item{5}{spot row number in a block/grid. Range [0,GROWS] in Z+.}
#   \item{6}{spot column number in a block/grid. Range [0,GCOLS] in Z+.}
#   \item{7}{spot (integrated) intensity in channel A corrected for
#     (median local) background. Note that saturated pixels are not used for calculations
#     of intensities as they can skew the results. Range [0,MAXAREA*65535]
#     in Z+.}
#   \item{8}{spot (integrated) intensity in channel B corrected for
#     (median local) background.  Note that saturated pixels are not used for calculations
#     of intensities as they can skew the results. Range [0,MAXAREA*65535]
#     in Z+.}
#
#  Column 9-17 exists only in full TAV files:\cr
#   \item{9}{spot mean ratio. Range [0,?] in R+. == (Col 7)/(Col 8), i.e. \bold{Redundant}.}
#
#   \item{10}{spot total area in pixels. Range [0,MAXAREA] in Z+.}
#   \item{11}{spot non-saturation factor. This measure shows the percentage 
#     of non-saturated pixels in the spot used for integration, i.e. values 
#     close to zero means a lot of saturation and values close to one means 
#     little saturation. Saturated pixels are not used for calculations of 
#     intensities as they can skew the results. Range [0,1] in R.}
#   \item{12}{spot median ratio = (median(signalA)-median(bkgA))/(median(signalB)-median(bkgB)). Range [0,MAXRATIO] in R.}
#   \item{13}{spot mode ratio = mode(signalA)-mode(nkgA))/(mode(signalB)-mode(bkgB)). Range [0,MAXRATIO] in R.}
#   \item{14}{total spot background in channel A =(median(BkgA)* spotarea). This background was subtracted when the spot intensity was calculated. The background is calculated by measuring the average local background and multiplying it by the spot area. Range [0,MAXAREA*65535] in Z+.}
#   \item{15}{total spot background in channel B =(median(BkgB)* spotarea). This background was subtracted when the spot intensity was calculated. The background is calculated by measuring the average local background and multiplying it by the spot area. Range [0,MAXAREA*65535] in Z+.}
#   \item{16}{spot flag in channel A. This flag is set by QC filter. Range: \{A,B,C,X,Y,Z\}.}
#   \item{17}{spot flag in channel B. This flag is set by QC filter.Range: \{A,B,C,X,Y,Z\}.}
#
#   The foreground estimates for channel A (B) can be obtained as the sum of column 7 and 14 (column 8 and 14).
# }
#
# \section{Information about spot QC flags}{
#   Flag values are generated based on next conditions: 
#    \item{A}{number of non-saturated pixels in spot is 0. 
#             \bold{Redundant}}
#    \item{B}{number of non-saturated pixels in spot is between 0 and 50.
#             \bold{Redundant}}
#    \item{C}{number of non-saturated pixels in spot is more then 50.
#             \bold{Redundant}}
#    \item{X}{spot was detected and rejected based on spot shape and spot
#             intensity relative to surrounding background.}
#    \item{Y}{spot background is higher than spot intensity. 
#             \bold{Redundant}}
#    \item{Z}{spot was not detected by the program. 
#             \bold{Redundant}}
# }
#
# @author
#
# \references{
#   TIGR Spotfinder Software, The Institute for Genomic Research, USA,
#   \url{http://www.tigr.org/software/tm4/}
# }
#
# \examples{\dontrun{
#   sf <- SpotfinderData$read("sf111.tav", path=system.file("data-ex", package="aroma"))
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
#   rg <- getRGData(ma)
#   plot(rg)
#   lowessCurve(rg, lwd=2, gridwise=TRUE)
#
#   # Plot M vs A with a lowess line through the data points
#   plot(ma)
#   lowessCurve(ma, lwd=2, gridwise=TRUE)
#
#   # Plot spatial
#   plotSpatial(ma)
# }}
#*/#########################################################################
setConstructorS3("SpotfinderData", function(layout=NULL) {
  extend(MicroarrayData(layout=layout), "SpotfinderData", 
    log.base = list()
  )
})

setMethodS3("append", "SpotfinderData", function(this, other) {
  NextMethod("append");
  invisible(this);
})



#########################################################################/**
# @RdocMethod readOneFile
#
# @title "Reads TAV files generated by Spotfinder"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filenames}{The filename(s) of the TAV/Spotfinder file(s) to be read.}
#   \item{path}{The path(s) to the file(s).}
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
#  Returns a @see "SpotfinderData" object containing one or several slides.
# }
#
# @author
#
# \examples{\dontrun{# For an example see help(SpotfinderData).}}
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("readOneFile", "SpotfinderData", function(this, filename, path=NULL, sep="auto", skip="auto", verbose=TRUE) {
  filename <- Arguments$getReadablePathname(filename, path);  

  if (verbose) cat("Reading file ", filename, "...", sep="");

  # TAV files are always tab delimeted headerless ASCII files.
  if (skip == "auto" || sep == "auto") {
    skip <- 0;
    sep <- "\t";
  }
  
  # Read the data with the header
  # The header/data *might* be quoted with " or ', but could be unquoted too.
  df <- read.table(file=filename, header=FALSE, sep=sep, quote="\"'", skip=skip, strip.white=TRUE, blank.lines.skip=TRUE);

  # Get version/type: 
  # "small" - "The original standard version of TAV file has 8 columns."
  # "full"  - Full version of TAV file is generated by Export ALL to TAV."
  if (ncol(df) == 8) {
    version <- "original";
  } else if (ncol(df) == 17) {
    version <- "full";
  } else {
    warning("Unknown number of columns read. It should be 8 (original TAV) or 17 (full TAV).");
  }
  
  # Rename the field names to standard field names.
  # Column 1-8 are always there
  colnames(df)[1:8] <- c("slideRow", "slideColumn", "metaRow", "metaColumn", "spotRow", "spotColumn", "integralA", "integralB");
  if (version == "full") {
    colnames(df)[9:17] <- c("meanRatio", "spotArea", "saturationFactor", "medianRatio", "modeRatio", "bkgA", "bkgB", "flagA", "flagB");
  }
  
  # Get the layout.
  gridfields <- c("metaRow", "metaColumn", "spotRow", "spotColumn");
  layout <- SpotfinderData$extractLayout(df[gridfields]);

  # Create 
  spot <- SpotfinderData(layout=layout)
  spot$.fieldNames <- names(df);
  for (field in names(df)) {
    spot[[field]] <- as.matrix(df[[field]]);
    df[[field]] <- NULL; # Save memory
  }

  if (verbose) cat("ok\n", sep="");
  
  rm(df); gc(); # To minimize memory usage!
  
  spot;
}, private=TRUE, static=TRUE);



#########################################################################/**
# @RdocMethod write
#
# @title "Write a SpotfinderData object to file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file to be written.}
#   \item{path}{The path to the file.}
#   \item{slide}{Slide to be saved.}
# }
#
# \value{Returns nothing.}
#
# @author
#
# \examples{\dontrun{
#   sf <- SpotfinderData$read("sf111.tav", path=system.file("data-ex", package="aroma"))
#
#   # Write the SpotfinderData object to a temporary file.
#   filename <- paste(tempfile("SpotfinderData"), ".tav", sep="")
#   write(sf, filename)
#
#   sf2 <- SpotfinderData$read(filename)
#   print(equals(sf, sf2))  # TRUE
#
#   unlink(filename)
# }}
#
# \seealso{
#   To read one or more SpotfinderData files at once
#   see @seemethod "read".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("write", "SpotfinderData", function(this, filename, path=NULL, slide=NULL, overwrite=FALSE, verbose=FALSE) {
  filename <- Arguments$getWritablePathname(filename, path, mustNotExist=!overwrite);  

  slide <- validateArgumentSlide(this, slide=slide);

  df <- extract(this, slides=slide);
  write.table(df, file=filename, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=FALSE);
});




#########################################################################/**
# @RdocMethod read
#
# @title "Reads several Spotfinder/TAV files into a SpotfinderData object"
#
# \description{
#  @get "title".
#  The acronym TAV stands for \emph{TIGR Array Viewer}. 
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
# \value{Returns a @see "SpotfinderData" object.}
#
# @author
#
# \examples{\dontrun{# For an example see help(SpotfinderData).}}
#
# \seealso{
#   For pattern formats see @see "base::list.files".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("read", "SpotfinderData", function(this, filename=NULL, path=NULL, pattern=NULL, verbose=TRUE, ...) {
  if (is.null(filename) && is.null(pattern))
    throw("Either 'filename' or 'pattern' must be specified.");
  if (!is.null(filename) && !is.null(pattern))
    throw("Only one of 'filename' and 'pattern' can be specified.");

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
    tmp <- SpotfinderData$readOneFile(filename=filename[k], path=path[k], verbose=verbose, ...);
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



setMethodS3("extractLayout", "SpotfinderData", function(this, griddata) {
  metaRow <- max(griddata[,1]);
  metaCol <- max(griddata[,2]);
  spotRow <- max(griddata[,3]);
  spotCol <- max(griddata[,4]);
  Layout(metaRow, metaCol, spotRow, spotCol);
}, static=TRUE);





#########################################################################/**
# @RdocMethod getRawData
#
# @title "Gets the raw intensites from the SpotfinderData structure"
#
# \description{
#  Extracts the red and green spot intensitites (both foreground and background)
#  from the SpotfinderData object and returns a @see "RawData" object.
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{Specifying which slides to be extracted. If @NULL, 
#     all slides are considered.}
# }
#
# \value{
#   Returns a @see "RawData" object containing the specified slides.
#   Note that the returned foreground and background signals for 
#   Spotfinder data is \emph{not} necessarily @integers, 
#   e.g. \code{R==0.32} is possible.
# }
#
# \details{
#   The R and Rb channels will come from the *A* fields, and
#   the G and Gb channels will come from the *B* fields.
#   To swap the channels use dyeSwap().
# }
#
# \examples{\dontrun{# For an example see help(SpotfinderData).}}
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getRawData", "SpotfinderData", function(this, slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  fieldNames <- getFieldNames(this);
  
  area <- this[["spotArea"]][,slides];
  Rb   <- this[["bkgA"]][,slides] / area;
  Gb   <- this[["bkgB"]][,slides] / area;
  R    <- (this[["integralA"]][,slides] + this[["bkgA"]][,slides]) / area;
  G    <- (this[["integralB"]][,slides] + this[["bkgB"]][,slides]) / area;
  rm(area);

  RawData(R=R, G=G, Rb=Rb, Gb=Gb, layout=getLayout(this), 
                                        extras=this$.extras)
})

setMethodS3("as.RawData", "SpotfinderData", function(this, ...) {
  getRawData(this, ...);
})


setMethodS3("getBackground", "SpotfinderData", function(this, which=c("mean")) {
  which <- match.arg(which);

  if (which == "mean") {
    fields <- c("bkgA", "bkgB")
  }

  # Assert that the fields do really exist.
  if (!all(fields %in% getFields(this))) {
    throw("The background estimates ", paste(fields, collapse=" and "), 
          " is not part of this ", data.class(this), " object.");
  }

  area <- this[["spotArea"]];
  RGData(R=this[[fields[1]]]/area, G=this[[fields[2]]]/area, layout=getLayout(this));
})



setMethodS3("getForeground", "SpotfinderData", function(this, which=c("mean")) {
  which <- match.arg(which);

  if (which == "mean") {
    fields <- c("integralA", "integralB")
  }

  # Assert that the fields do really exist.
  if (!all(fields %in% getFields(this))) {
    throw("The foreground estimates ", paste(fields, collapse=" and "), 
          " is not part of this ", data.class(this), " object.");
  }

  area <- this[["spotArea"]];
  R    <- (this[[fields[1]]] + this[["bkgA"]]) / area;
  G    <- (this[[fields[2]]] + this[["bkgB"]]) / area;
  rm(area);

  RGData(R=R, G=G, layout=getLayout(this));
})




#########################################################################/**
# @RdocMethod plotSpatial
#
# @title "Creates a spatial plot of a slide"
#
# \description{
#  @get "title".
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
# \examples{\dontrun{# For an example see help(SpotfinderData).}}
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("plotSpatial", "SpotfinderData", function(this, what="meanRatio", slide=1, pch=20, yaxt=NULL, ...) {
  NextMethod("plotSpatial", this, what=what, slide=slide, pch=pch, yaxt=NULL, ...);
})

setMethodS3("plotSpatial3d", "SpotfinderData", function(this, field="meanRatio", ...) {
  NextMethod("plotSpatial3d", this, field=field, ...);
})





setMethodS3("normalizeGenewise", "SpotfinderData", function(this, fields=NULL, bias=0, scale=1, ...) {
  if (is.null(fields)) {
    fields <- getFieldNames(this);
    exclude <- c("slideRow", "slideColumn", "metaRow", "metaColumn", 
                 "spotRow", "spotColumn", "flagA", "flagB");
    fields <- setdiff(fields, exclude);
  }
  NextMethod("normalizeGenewise", this, fields=fields, bias=bias, scale=scale, ...);
})


setMethodS3("getArea", "SpotfinderData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  as.matrix(this[["spotArea"]][include,slides]);
})


setMethodS3("getCircularity", "SpotfinderData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  throw("Sorry, but data from Spotfinder does not contain information such that the shape/circularity/perimiter can be calculated.");
})


setMethodS3("getSNR", "SpotfinderData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  throw("Sorry, but data from Spotfinder does not contain information such that the signal-to-noise ratio can be calculated.");
})


setMethodS3("getPerimeter", "SpotfinderData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  throw("Sorry, but data from Spotfinder does not contain information such that the shape/circularity/perimiter can be calculated.");
})


setMethodS3("getBgArea", "SpotfinderData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  throw("Sorry, but data from Spotfinder does not contain information about the number of pixels used in the background estimation.");
})



############################################################################
# HISTORY:
# 2005-10-21
# o Replace 'overwrite' arguments with 'mustNotExist' in calls to Arguments. 
# 2005-07-19
# o Replaced all path="" arguments to path=NULL.
# 2005-06-11
# o Making use of Arguments in R.utils.
# 2005-05-03
# o Updated regular expressions.
# 2004-08-15
# o Renamed as.RawData() to getRawData().
# 2004-03-09
# o BUG FIX: write() was mistakenly accepting more than once slide at the
#   time. It was a typo and now an error is thrown is more slides are given.
# 2003-05-14
# o BUG FIX: getForeground() and as.RawData() did by mistake return fg-2*bg
#   instead of fg as the foreground signal.
# 2003-04-12
# o Added getBackground() and getForeground().
# 2003-02-25
# o read(), write(), as.RawData(), plotSpatial(), plotSpatial3d(), extract()
#   seems to work.
# o Created from SpotData.
############################################################################
