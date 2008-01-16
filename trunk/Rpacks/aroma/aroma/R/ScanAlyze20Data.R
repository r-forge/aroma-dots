#########################################################################/**
# @RdocClass ScanAlyze20Data
#
# @title "The ScanAlyze20Data (ScanAlyze v2.0) class"
#
# @synopsis
#
# \description{
#  @classhierarchy
#
#  Creates an empty ScanAlyze20Data object.
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
# @author
#
# \references{
#   [1] \url{http://array.mbb.yale.edu/analysis/manual/inputfiles/}
# }
#*/#########################################################################
setConstructorS3("ScanAlyze20Data", function(layout=NULL) {
  extend(ScanAlyzeData(layout=layout), "ScanAlyze20Data")
})




#########################################################################/**
# @RdocMethod write
#
# @title "Writes a ScanAlyze v2.0 Results Data file (NOT IMPLEMENTED)"
#
# \description{
#  @get "title". 
# 
#  Not implemented yet!
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file to be read.}
#   \item{path}{The path to the file.}
#   \item{slides}{An @integer specifying which slides to be written to file. Currently, only one slide at the time can be written.}
#   \item{...}{Arguments passed to \code{write.table}.}
# }
#
# \value{Returns nothing.}
#
# @author
#
# \seealso{
#   To read one or more ScanAlyze v2.0 Results files see @seemethod "read".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("write", "ScanAlyze20Data", function(this, filename, path=NULL, slides=NULL, row.names=FALSE, ..., verbose=FALSE) {
  throw("Sorry, but writing ScanAlyze v2.0 is still not supported.");
});



setMethodS3("getKnownHeadersWithType", "ScanAlyze20Data", function(static) {
  # This methods returns a list of all known headers and their known
  # R datatype to be used when read by read.table() etc. 

  # NOTE: Some datatypes are currently unknown and are therefore set to
  # the read generic type NA, which is always understood. Hence, for
  # writing some of these types has to be updated before writing.

  # Regular expressions of headers. Watch out for "."'s and "+"'s!

  # ScanAlyze v 2.0 
  knownHeaders <- c(
    "SPOT"            ="integer", 
    "NAME"            ="character", 
    "Clone ID"        ="character", 
    "Gene Symbol"     ="character", 
    "Gene Name"       ="character", 
    "Cluster ID"      ="character", 
    "Accession"       ="character", 
    "Preferred name"  ="character", 
    "Locuslink ID"    ="character", 
    "SUID"            =NA, 
    "CH1I_MEAN"       =NA, 
    "CH1D_MEDIAN"     =NA, 
    "CH1I_MEDIAN"     =NA, 
    "CH1_PER_SAT"     =NA, 
    "CH1I_SD"         =NA, 
    "CH1B_MEAN"       =NA, 
    "CH1B_MEDIAN"     =NA, 
    "CH1B_SD"         =NA, 
    "CH1D_MEAN"       =NA, 
    "CH2I_MEAN"       =NA, 
    "CH2D_MEAN"       =NA, 
    "CH2D_MEDIAN"     =NA, 
    "CH2I_MEDIAN"     =NA, 
    "CH2_PER_SAT"     =NA, 
    "CH2I_SD"         =NA, 
    "CH2B_MEAN"       =NA, 
    "CH2B_MEDIAN"     =NA, 
    "CH2B_SD"         =NA, 
    "CH2BN_MEDIAN"    =NA, 
    "CH2DN_MEAN"      =NA, 
    "CH2IN_MEAN"      =NA, 
    "CH2DN_MEDIAN"    =NA, 
    "CH2IN_MEDIAN"    =NA, 
    "CORR"            =NA, 
    "DIAMETER"        =NA, 
    "FLAG"            =NA, 
    "LOG_RAT2N_MEAN"  =NA, 
    "LOG_RAT2N_MEDIAN"=NA, 
    "PIX_RAT2_MEAN"   =NA, 
    "PIX_RAT2_MEDIAN" =NA, 
    "PERGTBCH1I_1SD"  =NA, 
    "PERGTBCH1I_2SD"  =NA, 
    "PERGTBCH2I_1SD"  =NA, 
    "PERGTBCH2I_2SD"  =NA, 
    "RAT1_MEAN"       =NA, 
    "RAT1N_MEAN"      =NA, 
    "RAT2_MEAN"       =NA, 
    "RAT2_MEDIAN"     =NA, 
    "RAT2_SD"         =NA, 
    "RAT2N_MEAN"      =NA, 
    "RAT2N_MEDIAN"    =NA, 
    "REGR"            =NA, 
    "SUM_MEAN"        =NA, 
    "SUM_MEDIAN"      =NA, 
    "TOT_BPIX"        =NA, 
    "TOT_SPIX"        =NA, 
    "X_COORD"         =NA, 
    "Y_COORD"         =NA, 
    "TOP"             ="integer", 
    "BOT"             ="integer", 
    "LEFT"            ="integer", 
    "RIGHT"           ="integer", 
    "SECTOR"          ="integer", 
    "SECTORROW"       ="integer", 
    "SECTORCOL"       ="integer", 
    "SOURCE"          =NA, 
    "PLATEID"         =NA, 
    "PROW"            =NA, 
    "PCOL"            =NA, 
    "FAILED"          =NA, 
    "IS_VERIFIED"     =NA, 
    "IS_CONTAMINATED" =NA, 
    "LUID"            =NA
  );

  # ScanAlyze v?.?
  knownHeaders <- c(knownHeaders, 
    "GENE NAME"="character",
    "CHROMOSOME"="character",
    "STRAND"="character",
    "Beginning Coordinate"=NA, 
    "Ending Coordinate"=NA, 
    "MOLECULAR_FUNCTION"=NA,
    "BIOLOGICAL_PROCESS"=NA,
    "DS_GENE_1"=NA, 
    "DS_GENE_2"=NA, 
    "Forw Biol Primer Seq"=NA, 
    "Rev Biol Primer Seq"=NA
  )

  knownHeaders;
}, private=TRUE, static=TRUE);




##############################################################################
#  ScanAlyze 2.0 file format
#  
#  The number of descriptor columns (SPOT, NAME, etc) are variable.  
#  The data columns start wtih the heading "CH1I_MEAN".
#  
#  Each column must be separated by tab
#  ---------------------------------------------
#  
#  col datatype
#  0-x: variable number of descriptor columns
#  	 field           # what it means
#  + 0: CH1I_MEAN        # Uncorrected mean pixel intensity
#  + 1: CH1D_MEDIAN      # Median background subtracted intensities
#  + 2: CH1I_MEDIAN      # CH1 Intensity median
#  + 3: CH1_PER_SAT      # Percent saturation in CH1 of pixels
#  + 4: CH1I_SD          # spot intensity standard deviation
#  + 5: CH1B_MEAN        # Background mean intensity
#  + 6: CH1B_MEDIAN      # Background median intensity
#  + 7: CH1B_SD          # Background intensity standard deviation
#  + 8: CH1D_MEAN        # Ch1 Net intensity mean

#  + 9: CH2I_MEAN        # Uncorrected mean pixel intensity
#  +11: CH2D_MEDIAN      # Ch2 Net intensity median
#  +12: CH2I_MEDIAN      # Ch2 Intensity median
#  +13: CH2_PER_SAT      # Ch2 percent saturation of pixels
#  +14: CH2I_SD          # standard deviation of ch2
#  +15: CH2B_MEAN        # Background mean intensity
#  +16: CH2B_MEDIAN      # Background median intensity
#  +17: CH2B_SD          # Background intensity standard deviation
#  +10: CH2D_MEAN        # Ch2 Net Intensity Mean

#  +18: CH2BN_MEDIAN     # Ch2 normalized background median
#  +19: CH2DN_MEAN       # ch2 normalized net mean
#  +20: CH2IN_MEAN       # ch2 normalized intensity mean
#  +21: CH2DN_MEDIAN     # normalized ch2 net median
#  +22: CH2IN_MEDIAN     # normalized ch2 intensity median

#  +23: CORR             # Correlation between channel1 and channel2
#                        # pixels within spot
#  +24: DIAMETER         # spot diameter
#  +25: FLAG             # User defined spot flag
#  +26: LOG_RAT2N_MEAN   # log base 2 R/G normalized ratio mean
#  +27: LOG_RAT2N_MEDIAN # log base 2 R/G normalized ratio median
#  +28: PIX_RAT2_MEAN    # R/G mean per pixel
#  +29: PIX_RAT2_MEDIAN  # R/G median per pixel
#  +30: PERGTBCH1I_1SD   # % Ch1 pixels > BG + 1SD
#  +31: PERGTBCH1I_2SD   # % Ch1 pixels > BG + 2SD
#  +32: PERGTBCH2I_1SD   # % Ch2 pixels > BG + 1SD
#  +33: PERGTBCH2I_2SD   # % Ch2 pixels > BG + 2SD
#  +34: RAT1_MEAN        # G/R mean
#  +35: RAT1N_MEAN       # G/R normalized mean
#  +36: RAT2_MEAN        # R/G mean
#  +37: RAT2_MEDIAN      # R/G median
#  +38: RAT2_SD          # std dev of pixel intensity ratios
#  +39: RAT2N_MEAN       # R/G normalized mean
#  +40: RAT2N_MEDIAN     # R/G normalized median
#  +41: REGR             # Slope of linear regression line to data
#  +42: SUM_MEAN         # sum of mean intensities
#  +43: SUM_MEDIAN       # sum of median intensities
#  +44: TOT_BPIX         # Number of pixels in background
#  +45: TOT_SPIX         # Number of pixels in spot
#  +46: X_COORD          # of spot in image
#  +47: Y_COORD          # of spot in image
#  +48: TOP              # Coordinates of box containing spot ellipse
#  +49: BOT              # Coordinates of box containing spot ellipse
#  +50: LEFT             # Coordinates of box containing spot ellipse
#  +51: RIGHT            # Coordinates of box containing spot ellipse
#  +52: SECTOR           # Grid number (left to right, top to bottom)
#  +53: SECTORROW        # Row of spot in grid
#  +54: SECTORCOL        # column of spot in grid
#  +55: SOURCE     
#  +56: PLATEID     
#  +57: PROW     
#  +58: PCOL     
#  +59: FAILED     
#  +60: IS_VERIFIED     
#  +61: IS_CONTAMINATED     
#  +62: LUID
#
# Source: http://array.mbb.yale.edu/analysis/manual/inputfiles/
##############################################################################
setMethodS3("readInternal", "ScanAlyze20Data", function(this, filename) {
  # This method assumes that there is a bunch of ! lines (comments), followed by 
  # a HEADER line, (that begins with SPOT). All following lines are then assumed 
  # to be Spot data lines (beginning with a number and then data entries 
  # according to the header) /Toby Dylan Hocking, 2003-10-07
  fh <- file(filename, "r");
  on.exit({ 
    if (!is.null(fh)) {
      close(fh);
      fh <- NULL;
    }
  });

  # Read until the header line (SPOT) is found.
  remark <- NULL;
  kk <- 0;
  isComment <- TRUE;
  while (isComment) {
    kk <- kk + 1;
    line <- scan(file=fh, what=character(0), sep="\t", quote="", quiet=TRUE, strip.white=TRUE, blank.lines.skip=TRUE, nlines=1);
    lineType <- line[1];
    isComment <- (regexpr("^!",lineType) != -1);
    isComment <- isComment || (regexpr("^\"!",lineType) != -1);
    if (isComment) {
      line <- gsub("^!", "", line);
      line <- gsub("=", "\t", line);
      remark <- c(remark, line);
    }
  }
  close(fh);
  fh <- NULL;

  if (lineType != "SPOT") {
    throw("ScanAlyze file format error: Expected a header line (starting with \"SPOT\") as the first occurance after comments.");
  } 

  nbrOfHeaderLines <- kk;

  header <- line;
  # Remove trailing blanks
  header <- gsub("[ ]*$", "", header);
  # Remove quotation marks
  header <- gsub("\"", "", header);

  knownHeaders <- getKnownHeadersWithType(this);

  # The below code is inspired by the GenePixData code.
  saVersion <- NA;
  unknownHeaders <- c();
  colClasses <- c();
  for (kk in seq(along=header)) {
    colNames <- names(knownHeaders);
    typeIdx <- which(sapply(seq(colNames), FUN=function(i) regexpr(colNames[i], header[kk]) == 1));
    if (length(typeIdx) == 0) {
      unknownHeaders <- c(unknownHeaders, kk);
      colClasses <- c(colClasses, NA);
    } else {
      typeIdx <- typeIdx[1];  # Bullet proof!
      colClasses <- c(colClasses, knownHeaders[typeIdx]);
    }
  }

  if (length(unknownHeaders) > 0) {
    unknown <- paste(header[unknownHeaders], collapse=", ");
    warning(paste(sep="", "ScanAlyze file format warning: Some of the field names are not recognized. It might be because it is a new/old version. The unknown fields are: ", unknown, "."));
  }

  fh <- file(filename, "r");
  on.exit({ 
    if (!is.null(fh)) {
      close(fh);
      fh <- NULL;
    }
  });
  nbrOfColumns <- length(header);
#  df <- scan(file=fh, what=rep(list(""), nbrOfColumns+1), sep="\t", quote="", quiet=TRUE, strip.white=TRUE, blank.lines.skip=TRUE, skip=nbrOfHeaderLines);
  close(fh);
  fh <- NULL;

  df <- read.table(filename, sep="\t", quote="", comment.char="", skip=nbrOfHeaderLines, colClasses=colClasses, header=FALSE, na=c("NA", "NaN", "Error"));
  names(df) <- header;

  # Convert all columns.
  for (kk in 1:nbrOfColumns) {
    x <- type.convert(as.character(df[[kk]]), na.strings=c("Inf","NA"));
#    if (!all(is.na(x)))
#      x[df[[kk]] == "Inf"] <- Inf;
    df[[kk]] <- x;
  }

   # Add a 'slide' column
   slide <- rep(1, length(df[[1]]));
   df <- c(list(slide), df);

  # Put names on columns
  names(df) <- c(header);
  names(df) <- c("slide", header);
  
  list(header=header, remark=remark, spot=as.data.frame(df), version="2.0");
}, private=TRUE, static=TRUE);




#########################################################################/**
# @RdocMethod read
#
# @title "Reads one or several ScanAlyze v2.0 files"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Same arguments as read() in ScanAlyzeData class.}
# }
#
# \value{Returns a @see "ScanAlyze20Data" object.}
#
# @author
#
# \references{
#   \url{http://array.mbb.yale.edu/analysis/manual/inputfiles/}
# }
#
# \seealso{
#   For pattern formats see @see "base::list.files".
#   @see "ScanAlyzeData".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("read", "ScanAlyze20Data", function(static, ...) {
  # Entering here 'static' is a ScanAlyze20Data object, which is then
  # explicitly passed to read() in the super class where it is used
  # to identify the correct readOneFile() method.
  read.ScanAlyzeData(static, ...);
}, static=TRUE);





############################################################################
#
#   S I G N A L    E S T I M A T E S
#
############################################################################
setMethodS3("getForeground", "ScanAlyze20Data", function(this, which=c("median", "mean", "auto"), slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  which <- match.arg(which);
  if (which == "auto")
    which <- "median";

  if (which == "mean") {
    fields <- c("CH2I.MEAN", "CH1I.MEAN");
  } else if (which == "median") {
    fields <- c("CH2I.MEDIAN", "CH1I.MEDIAN");
  }

  # Assert that the fields do really exist.
  if (!all(fields %in% getFields(this))) {
    throw("The foreground estimates ", paste(fields, collapse=" and "), 
          " is not part of this ", data.class(this), " object.");
  }

  rg <- RGData(R=this[[fields[1]]][,slides], G=this[[fields[2]]][,slides], layout=getLayout(this));
  if (all(is.na(rg$R) & is.na(rg$G))) {
    warning("Could not find any median foreground estimates. Revert to the mean.");
    getForeground(this, which="mean", slides=slides);
  } else {
    rg;
  }
})

setMethodS3("getForegroundSD", "ScanAlyze20Data", function(this, slides=NULL, index=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(index))
    index <- 1:nbrOfSpots(this);
  
  list(RSD=this[["CH2I.SD"]][index,slides], GSD=this[["CH1I.SD"]][index,slides]);
})


setMethodS3("getBackground", "ScanAlyze20Data", function(this, which=c("median", "mean", "auto"), slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  which <- match.arg(which);
  if (which == "auto")
    which <- "median";

  if (which == "mean") {
    fields <- c("CH2B.MEAN", "CH1B.MEAN");
  } else if (which == "median") {
    fields <- c("CH2B.MEDIAN", "CH1B.MEDIAN");
  }

  # Assert that the fields do really exist.
  if (!all(fields %in% getFields(this))) {
    throw("The background estimates ", paste(fields, collapse=" and "), 
          " is not part of this ", data.class(this), " object.");
  }

  RGData(R=this[[fields[1]]][,slides], G=this[[fields[2]]][,slides], layout=getLayout(this));
})


setMethodS3("getBackgroundSD", "ScanAlyze20Data", function(this, slides=NULL, index=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(index))
    index <- 1:nbrOfSpots(this);
  
  list(RSD=this[["CH2B.SD"]][index,slides], GSD=this[["CH1B.SD"]][index,slides]);
})



############################################################################
#
#   S P O T   G E O M E T R Y
#
############################################################################
setMethodS3("getArea", "ScanAlyze20Data", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  # "+45: TOT.SPIX         # Number of pixels in spot"
  as.matrix(this[["TOT.SPIX"]][include,slides]);
})


setMethodS3("getBgArea", "ScanAlyze20Data", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  # "+44: TOT.BPIX         # Number of pixels in background"
  as.matrix(this[["TOT.BPIX"]][include,slides]);
})


setMethodS3("getDiameter", "ScanAlyze20Data", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  this[["DIAMETER"]][include,slides];
})



setMethodS3("getGrid", "ScanAlyze20Data", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  this[["SECTOR"]][include,slides];
})

setMethodS3("getSpotRow", "ScanAlyze20Data", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  this[["SECTORROW"]][include,slides];
})

setMethodS3("getSpotColumn", "ScanAlyze20Data", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  this[["SECTORCOL"]][include,slides];
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
setMethodS3("getSpotPosition", "ScanAlyze20Data", function(this, slides=NULL, index=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (is.null(index)) {
    index <- 1:nbrOfSpots(this);
  } else if (any(index < 1) || any(index > nbrOfSpots(this))) {
    throw("Argument 'index' is out of range.");
  }
  
  xField <- "X.COORD";
  yField <- "Y.COORD";
  if (hasField(this, xField) && hasField(this, yField) &&
      !all(is.na(this[[xField]])) && !all(is.na(this[[yField]]))) {
    x <- this[[xField]][index,slides];
    y <- this[[yField]][index,slides];
    SpotPosition(x=x, y=y);
  } else {
    NextMethod("getSpotPosition");
  }
})




############################################################################
# HISTORY:
# 2005-07-19
# o Replaced all path="" arguments to path=NULL.
# 2003-10-27
# o Added support to read ScanAlyze v2.0 files. Still have to document it.
# o Since ScanAlyze v2.0 is soo different from before, I decided to make it
#   its own class ScanAlyze20Data.
#   Right now it inherits from ScanAlyzeData, because then one can just
#   override broken methods and reuse the others (=most).
# o Created from ScanAlyzeData.R.
#
############################################################################
