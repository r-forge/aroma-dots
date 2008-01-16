#########################################################################/**
# @RdocClass ScanAlyzeData
#
# @title "The ScanAlyzeData class"
#
# @synopsis
#
# \description{
#  @classhierarchy
#
#  Creates an empty ScanAlyzeData object.
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
#  A ScanAlyze file contains spot information for each spot on a single
#  microarray slide. It consists of a header followed by a unspecified 
#  number of rows. The header contains a number of field labels. Each row
#  starting with \code{SPOT} corresponds to one spot. Other rows starts 
#  with \code{HEADER} or \code{REMARK}. The spot fields are:\cr
#
#   \item{SPOT}{Unique index of spot in file. "Counting starts with 
#     grid 1, moves along row 1 from column 1 until the last column, then
#     advances to the next row; after all rows in grid 1 are assigned an 
#     index, counting proceeds to grid 2, etc.". Comment: This is the same
#     way as GenePix and Spot indices the spots.}
#   \item{GRID}{Grid number.}
#   \item{TOP, LEFT, BOT, RIGHT}{Coordinates of the box containing spot 
#     ellipse, in image coordinates.}
#   \item{ROW}{Row within grid.}
#   \item{COL}{Column within grid.}
#   \item{CH1I, CH2I}{Channel 1,2 - uncorrected mean of the foreground pixel
#     intensites.}
#   \item{CH1B, CH2B}{Channel 1,2 - median of the background pixel intensites.}
#   \item{CH1BA, CH2BA}{Channel 1,2 - mean of the background pixel intensites.}
#   \item{SPIX}{Number of pixels in the foreground, i.e. in the spot.}
#   \item{BGPIX}{Number of pixels in the background.}
#   \item{EDGE}{?}
#   \item{RAT2}{Median of (CH2PI-CH2B)/(CH1PI-CH1B) where CH1PI and CH2PI are
#     the single pixel intensities.}
#   \item{MRAT}{?}
#   \item{REGR}{?}  
#   \item{CORR}{The correlation between channel1 and channel2 pixels within
#     the spot.}
#   \item{LFRAT}{?}
#   \item{CH1GTB1, CH2GTB1}{Fraction of pixels in the spot greater than the
#     background (CH1B or CH2B).}
#   \item{CH1GTB2, CH2GTB2}{Fraction of pixels in the spot greater than 1.5 times the background (CH1B or CH2B).}
#   \item{CH1EDGEA, CH2EDGEA}{mean magnitude of the horizontal and vertical Sobel edge vectors within spot 1 and spot 2, respectively}
#   \item{FLAG}{User defined spot flag (default 0).}
#   \item{CH1KSD, CH2KSD}{The value of the Komogorov-Smirnov statistic [3] that assesses the likelihood that the spot pixel intensity distribution is drawn from the background distribution.}
#   \item{CH1KSP, CH2KSP}{The actual probabilities of the above statistic.}
# }
#
#
# @author
#
# \references{
#  [1] ScanAlyze Software, Eisen Lab, Lawrence Berkeley National Lab, 2001,
#   \url{http://rana.lbl.gov/EisenSoftware.htm}.\cr
#  [2] ScanAlyze User Manual, Michael Eisen, Stanford University, 1999,
#   \url{http://rana.lbl.gov/manuals/ScanAlyzeDoc.pdf}.\cr
#  [3] 14.3 Are Two Distributions Different? a sample chapter from 
#      Numerical Recipes in C: The Art of Scientific Computing,
#      1992, Cambridge University Press
#   \url{http://www.ulib.org/webRoot/Books/Numerical_Recipes/bookcpdf/c14-3.pdf}.\cr
# }
#
# \examples{
#   sa <- ScanAlyzeData$read("group4.dat", path=system.file("data-ex", package="aroma"))
#
#   # Get the foreground and the background (and the layout)
#   raw <- getRawData(sa)
#
#   # The the background corrected data
#   ma <- getSignal(raw, bgSubtract=FALSE)
#
#   # Plot M vs A with a lowess line through the data points
#   plot(ma, slide=1)
#   lowessCurve(ma, lwd=2, gridwise=TRUE)
# }
#
# \note{
#   The laser beam with wavelength 635nm is the red laser, and the one with
#   wavelength 532nm is the green laser.
# }
#*/#########################################################################
setConstructorS3("ScanAlyzeData", function(layout=NULL) {
  extend(MicroarrayData(layout=layout), "ScanAlyzeData", 
    remark = list()
  )
})


setMethodS3("append", "ScanAlyzeData", function(this, other) {
  NextMethod("append");
  if (inherits(other, "ScanAlyzeData"))
    this$remark <- c(this$remark, other$remark)
  else
    this$remark <- c(this$remark, NA);
  invisible(this);
})



#########################################################################/**
# @RdocMethod read
#
# @title "Reads data from a ScanAlyze results data file"
#
# \description{
#  @get "title". 
#  This method will also read ScanAlyze v2.0 files. In such cases a
#  ScanAlyze20Data object is returned!
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the ScanAlyze file to be read.}
#   \item{path}{The path(s) to the ScanAlyze file.}
#   \item{safe}{If @TRUE even obscure ScanAlyze formats should be possible
#     to read, but this method is \emph{really} slow! However, most 
#     likely if no one has modified the data file by hand, the standard
#     reading method could be used, which is \emph{much} faster.}
# }
#
# \value{
#  Returns a @see "ScanAlyzeData" object containing one or several slides.
# }
#
# \author{@get "author". Code for reading ScanAlyze v2.0 files was provided 
#  by Toby Dylan Hocking at UC Berkeley.}
#
# \examples{
#   sa <- ScanAlyzeData$read("group4.dat", path=system.file("data-ex", package="aroma"))
# }
#
# \references{
#   The data file in the example above was contributed by
#   Michael Stadler, Bioinformatics Group, Swiss Institute for
#   Experimental Cancer Research.
# }
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("readOneFile", "ScanAlyzeData", function(static, filename, path=NULL, safe=FALSE, verbose=FALSE) {
  filename <- Arguments$getReadablePathname(filename, path);  

  if (verbose) cat("Reading file ", filename, "...", sep="");

  # Support gzip'ed files too.
  if (regexpr("[.]gz$", filename) != -1) {
    tmpname <- tempfile();
    n <- gunzip(filename, tmpname);
    filename <- tmpname;
    on.exit(file.remove(tmpname));
  }
  

  read <- constructor <- NULL;
  tryCatch({
    if (safe) {
      read <- ScanAlyzeData$readInternalSafe(filename)
    } else {
      read <- ScanAlyzeData$readInternal(filename);
    }
    constructor <- ScanAlyzeData;
  }, error=function(ex) {
    # Try ScanAlyze v2.0 if not found. Superclasses could always
    # try its subclasses, but not the otherway around.
    read <<- ScanAlyze20Data$readInternal(filename);
    constructor <<- ScanAlyze20Data;
  })

  header <- read$header;
  remark <- read$remark;
  df <- read$spot;
  version <- read$version;
  rm(read);

  if (is.element("GRID", header)) {
    gridField    <- "GRID";
    spotRowField <- "ROW";
    spotColField <- "COL";
  } else if (is.element("SECTOR", header)) {
    # Hmm, here we are actually cheating, because we use v2.0 fields
    # in here.
    gridField    <- "SECTOR";
    spotRowField <- "SECTORROW";
    spotColField <- "SECTORCOL";
  } else {
    throw("No grid information found in ScanAlyze data file: ", filename);
  }

  nbrOfGrids <- max(as.integer(df[,gridField]), na.rm=TRUE);
  nspot.r    <- max(as.integer(df[,spotRowField]), na.rm=TRUE);
  nspot.c    <- max(as.integer(df[,spotColField]), na.rm=TRUE);

  gridSize <- nspot.r*nspot.c;

  nmissing <- gridSize*nbrOfGrids-nrow(df);
  if (nmissing > 0) {
    if (verbose) {
      cat("Filling in missing observations, i.e. excluded spot rows: ", nmissing, "\n");
    }

    # Allocate a new data frame with the same fields, but with *all* rows
    df2 <- df[rep(1,gridSize*nbrOfGrids),];

    # Fill in missing rows
    shortGrids <- NULL;
    spotRowIdx <- rep(1:nspot.r, each=nspot.c);
    spotColIdx <- rep(1:nspot.c, length=length(spotRowIdx));
    spotRowColIdx <- c(spotRowIdx,spotColIdx);
    spotRowColStr <- paste(spotRowIdx, spotColIdx, sep=",");
    for (grid in 1:nbrOfGrids) {
      gridIdx  <- which(as.integer(df[,gridField]) == grid);
      offset <- (grid-1)*gridSize; 
      if (length(gridIdx) == gridSize) {
        idx <- offset + 1:gridSize;
        df2[idx,] <- df[gridIdx,];
      } else {
        hasSpotRowIdx <- as.integer(df[,c(spotRowField)]);
        hasSpotColIdx <- as.integer(df[,c(spotColField)]);
        hasSpotRowColStr <- paste(hasSpotRowIdx, hasSpotColIdx, sep=",");
        match <- match(spotRowColStr, hasSpotRowColStr);
        missing <- is.na(match);
        # Set all missing rows to NA
        idx <- offset + which(missing);
        df2[idx,] <- NA;
        df2[idx,"SPOT"] <- idx;
        df2[idx,gridField] <- grid;
        df2[idx,spotRowField] <- spotRowIdx;
        df2[idx,spotColField] <- spotColIdx;
        # ...and copy the data to all non-missing rows
        idx <- offset + which(!missing);
        df2[idx,] <- df[gridIdx,];
#        missingPairs <- spotRowColIdx[missing];
#        msg <- paste(spotRowColStr[missing], collapse="; ");
#        msg <- paste("Grid #", grid, ": ", msg, sep="");
#        print(msg);
      }
    }
    df <- df2;
    rm(df2);
#    throw(sprintf("The number of rows in the ScanAlyze file is not equal to the expected %d grids with %dx%d (=%d) spots: %d", nbrOfGrids, nspot.r, nspot.c, gridSize*nbrOfGrids, nrow(df)));
  }

  # Get the index of the first spot in each grid.
  idx <- gridSize*(0:(nbrOfGrids-1)) + 1;
  # Get the position of the first spot in each grid.
  x <- df[idx,"LEFT"];
  # Get the distance in the x-dimension to the next grid.
  dx <- c(diff(x), -1);
  dx <- -sign(sign(dx)-1);
  ngrid.r <- sum(dx, na.rm=TRUE);
  ngrid.c <- nbrOfGrids / ngrid.r;

  # Assert that out calculation are correct.
  if (round(ngrid.c) != ngrid.c)
    warning("Number of grid rows seems to be wrong! Please report this error to the author of the package.\n");

  layout <- Layout(ngrid.r, ngrid.c, nspot.r, nspot.c)
  sa <- constructor(layout=layout)
  sa$version <- version;
  sa$remark <- list(remark);
  sa$.fieldNames <- names(df);
  for (field in names(df)) {
    sa[[field]] <- as.matrix(df[[field]]);
    df[[field]] <- NULL; # Save memory
  }

  if (verbose) cat("ok\n", sep="");
  
  rm(df); gc(); # To minimize memory usage!

  sa;
}, private=TRUE, static=TRUE);



#########################################################################/**
# @RdocMethod write
#
# @title "Write a ScanAlyze Results Data file"
#
# \description{
#  Writes the ScanAlyzeData object to a file with the ScanAlyze file format.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file to be read.}
#   \item{path}{The path to the file.}
#   \item{slide}{An @integer specifying which slides to be written to file. Currently, only one slide at the time can be written.}
#   \item{...}{Arguments passed to \code{write.table}.}
# }
#
# \value{Returns nothing.}
#
# @author
#
# \examples{
#   sa <- ScanAlyzeData$read("group4.dat", path=system.file("data-ex", package="aroma"))
#
#   # Writes the ScanAlyzeData to a temporary file. Note that this
#   # file won't be exactly the same since a few lines, specifying
#   # for instance the creator of the file, will be added. The data,
#   # however, will be the same.
#   filename <- paste(tempfile("ScanAlyzeData"), ".dat", sep="")
#   write(sa, filename)
#
#   sa2 <- ScanAlyzeData$read(filename)
#   print(equals(sa, sa2))  # TRUE
#
#   unlink(filename)
# }
#
# \references{
#   The data file in the example above was contributed by
#   Michael Stadler, Bioinformatics Group, Swiss Institute for
#   Experimental Cancer Research.
# }
#
# \seealso{
#   To read one or more ScanAlyze Results files
#   see @seemethod "read".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("write", "ScanAlyzeData", function(this, filename, path=NULL, slide=NULL, overwrite=FALSE, row.names=FALSE, ..., verbose=FALSE) {
  writeVersion10 <- function(fh, slide=slide) {
    # Write the "HEADER" line.
    cat(file=fh, "HEADER\t", paste(getFieldNames(this), collapse="\t"), "\n", sep="");
    
    # Write the "REMARK" line(s).
    remark <- this$remark[[slide]];
    if (!is.element("CREATOR", names(remark)))
      remark <- c(remark, CREATOR=paste("com.braju.sma, http://www.braju.com/R/,", date()));
    
    for (k in seq(remark)) {
      cat(file=fh, "REMARK\t", names(remark)[k], "\t", remark[k], "\n", sep="");
    }
    
    df <- extract(this, slides=slide);
    # Adds the "SPOT" column.
    df <- cbind(rep("SPOT", nrow(df)), df);
    write.table(df, file=fh, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE, ...);
  } # writeVersion10()

  
  filename <- Arguments$getWritablePathname(filename, path, mustNotExist=!overwrite);  

  slide <- validateArgumentSlide(this, slide=slide);

  fh <- file(filename, "w");
  on.exit(close(fh));
  
  if (this$version == "1.0") {
    writeVersion10(fh, slide=slide);
  } else {
    throw("Sorry, but currently only support for writing ScanAlyze v1.0 to file: ", this$version);
  }
});



setMethodS3("getKnownHeadersWithType", "ScanAlyzeData", function(static) {
  # This methods returns a list of all known headers and their known
  # R datatype to be used when read by read.table() etc. 

  # NOTE: Some datatypes are currently unknown and are therefore set to
  # the read generic type NA, which is always understood. Hence, for
  # writing some of these types has to be updated before writing.

  # Regular expressions of headers. Watch out for "."'s and "+"'s!
  # Note that the "Log Ratio" field sometimes contains "Error", i.e. we
  # can't use "double".

  # ScanAlyze v 1?
  knownHeaders <- c(
    "SPOT"="integer", 
    "GRID"="integer",
    "TOP"="integer",
    "LEFT"="character",
    "BOT"="character",
    "RIGHT"="integer",
    "ROW"="integer",
    "COL"="integer",
    "CH1I"=NA,
    "CH1B"=NA,
    "CH1AB"=NA, 
    "CH2I"=NA, 
    "CH2B"=NA, 
    "CH2AB"=NA, 
    "SPIX"=NA, 
    "BGPIX"=NA, 
    "EDGE"=NA, 
    "RAT2"=NA, 
    "MRAT"=NA, 
    "REGR"=NA, 
    "CORR"=NA, 
    "LFRAT"="double", 
    "CH1GTB1"=NA, 
    "CH2GTB1"=NA, 
    "CH1GTB2"=NA, 
    "CH2GTB2"=NA, 
    "CH1EDGEA"=NA, 
    "CH2EDGEA"=NA, 
    "FLAG"=NA, 
    "CH1KSD"=NA, 
    "CH1KSP"=NA, 
    "CH2KSD"=NA, 
    "CH2KSP"=NA
  );

  # ScanAlyze v 2.0 
  knownHeaders <- c(knownHeaders, 
    "SPOT"="integer", 
    "NAME"="character", 
    "Clone ID"="character", 
    "Gene Symbol"="character", 
    "Gene Name"="character", 
    "Cluster ID"="character", 
    "Accession"="character", 
    "Preferred name"="character", 
    "Locuslink ID"="character", 
    "SUID"=NA, 
    "CH1I_MEAN"=NA, 
    "CH1D_MEDIAN"=NA, 
    "CH1I_MEDIAN"=NA, 
    "CH1_PER_SAT"=NA, 
    "CH1I_SD"=NA, 
    "CH1B_MEAN"=NA, 
    "CH1B_MEDIAN"=NA, 
    "CH1B_SD"=NA, 
    "CH1D_MEAN"=NA, 
    "CH2I_MEAN"=NA, 
    "CH2D_MEAN"=NA, 
    "CH2D_MEDIAN"=NA, 
    "CH2I_MEDIAN"=NA, 
    "CH2_PER_SAT"=NA, 
    "CH2I_SD"=NA, 
    "CH2B_MEAN"=NA, 
    "CH2B_MEDIAN"=NA, 
    "CH2B_SD"=NA, 
    "CH2BN_MEDIAN"=NA, 
    "CH2DN_MEAN"=NA, 
    "CH2IN_MEAN"=NA, 
    "CH2DN_MEDIAN"=NA, 
    "CH2IN_MEDIAN"=NA, 
    "CORR"=NA, 
    "DIAMETER"=NA, 
    "FLAG"=NA, 
    "LOG_RAT2N_MEAN"=NA, 
    "LOG_RAT2N_MEDIAN"=NA, 
    "PIX_RAT2_MEAN"=NA, 
    "PIX_RAT2_MEDIAN"=NA, 
    "PERGTBCH1I_1SD"=NA, 
    "PERGTBCH1I_2SD"=NA, 
    "PERGTBCH2I_1SD"=NA, 
    "PERGTBCH2I_2SD"=NA, 
    "RAT1_MEAN"=NA, 
    "RAT1N_MEAN"=NA, 
    "RAT2_MEAN"=NA, 
    "RAT2_MEDIAN"=NA, 
    "RAT2_SD"=NA, 
    "RAT2N_MEAN"=NA, 
    "RAT2N_MEDIAN"=NA, 
    "REGR"=NA, 
    "SUM_MEAN"=NA, 
    "SUM_MEDIAN"=NA, 
    "TOT_BPIX"=NA, 
    "TOT_SPIX"=NA, 
    "X_COORD"=NA, 
    "Y_COORD"=NA, 
    "TOP"="integer", 
    "BOT"="integer", 
    "LEFT"="integer", 
    "RIGHT"="integer", 
    "SECTOR"="integer", 
    "SECTORROW"="integer", 
    "SECTORCOL"="integer", 
    "SOURCE"=NA, 
    "PLATEID"=NA, 
    "PROW"=NA, 
    "PCOL"=NA, 
    "FAILED"=NA, 
    "IS_VERIFIED"=NA, 
    "IS_CONTAMINATED"=NA, 
    "LUID"=NA
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



setMethodS3("readInternalSafe", "ScanAlyzeData", function(this, filename) {
  headerV2.44 <- c("SPOT", "GRID", "TOP", "LEFT", "BOT", "RIGHT", "ROW", "COL", "CH1I", "CH1B", "CH1AB", "CH2I", "CH2B", "CH2AB", "SPIX", "BGPIX", "EDGE", "RAT2", "MRAT", "REGR", "CORR", "LFRAT", "CH1GTB1", "CH2GTB1", "CH1GTB2", "CH2GTB2", "CH1EDGEA", "CH2EDGEA", "FLAG", "CH1KSD", "CH1KSP", "CH2KSD", "CH2KSP");

  # This method assumes nothing about the file structure more than there is *one*
  # HEADER line, some REMARK lines and some SPOT line, in any order. It is safer,
  # but also *very* slow since it has to read line by line.
  fh <- file(filename, "r");
  
  header <- NULL;
  remark <- NULL;
  df <- NULL;
  isAtEndOfFile <- FALSE;
  while (!isAtEndOfFile) {
    line <- scan(file=fh, what=character(0), sep="\t", quote="", quiet=TRUE, strip.white=TRUE, blank.lines.skip=TRUE, nlines=1);
    lineType <- line[1];
    if (lineType == "SPOT") {
      data <- line[-1];
      data <- as.numeric(data);
      df <- rbind(df, data);
    } else if (lineType == "HEADER") {
      if (is.null(header)) {
        header <- line[-1];
        if (any(header != headerV2.44)) {
          warning("ScanAlyze file format warning: The HEADER of the ScanAlyze data file is not recognized. Please be aware of this!. Some of the field names are not recognized.");
        }
      } else {
        throw("ScanAlyze file format error: There are more than one row containing a HEADER.");
      }
    } else if (lineType == "REMARK") {
      remark <- c(remark, paste(line[-c(1,2)], sep="\t"));
      names(remark)[length(remark)] <- line[2];
    } else if (length(line) == 0) {
      isAtEndOfFile <- TRUE;
    } else {
      throw("ScanAlyze file format error: There is a line starting with \"", lineType, "\", which is an unknown line type.");
    }
  } # while(!isAtEndOfFile)
  close(fh);
  dimnames(df) <- list(NULL, header);
  list(header=header, remark=remark, spot=df, version="1.0");
}, private=TRUE, static=TRUE);



setMethodS3("readInternal", "ScanAlyzeData", function(this, filename) {
  headerV2.44 <- c("HEADER", "SPOT", "GRID", "TOP", "LEFT", "BOT", "RIGHT", "ROW", "COL", "CH1I", "CH1B", "CH1AB", "CH2I", "CH2B", "CH2AB", "SPIX", "BGPIX", "EDGE", "RAT2", "MRAT", "REGR", "CORR", "LFRAT", "CH1GTB1", "CH2GTB1", "CH1GTB2", "CH2GTB2", "CH1EDGEA", "CH2EDGEA", "FLAG", "CH1KSD", "CH1KSP", "CH2KSD", "CH2KSP");
  
  # This method assumes that there is a HEADER line, followed by some REMARK lines
  # and then *everything else* following afterwards are SPOT lines. This is probably
  # the normal case, but if someone puts a REMARK after the first SPOT line, this
  # method will give an error.
  fh <- file(filename, "r");

  line <- scan(file=fh, what=character(0), sep="\t", quote="", quiet=TRUE, strip.white=TRUE, blank.lines.skip=TRUE, nlines=1);
  lineType <- line[1];
  if (lineType != "HEADER") {
    throw("ScanAlyze file format error: Expected a HEADER line as the first occurance. Try to read the file with safe=TRUE.");
  } else if (all(line == headerV2.44)) {
    header <- line[-1];
  } else {
    header <- line[-1];
    warning("ScanAlyze file format warning: The HEADER of the ScanAlyze data file is not recognized. Please be aware of this!");
  }
  
  # Read until the first SPOT line is found.
  remark <- NULL;
  kk <- 0;
  isSpot <- FALSE;
  while (!isSpot) {
    line <- scan(file=fh, what=character(0), sep="\t", quote="", quiet=TRUE, strip.white=TRUE, blank.lines.skip=TRUE, nlines=1);
    if (length(line) > 0) {
      lineType <- line[1];
      if (lineType == "REMARK") {
  	remark <- c(remark, paste(line[-c(1,2)], sep="\t"));
  	names(remark)[length(remark)] <- line[2];
      }
      isSpot <- (lineType == "SPOT");
    }
    kk <- kk + 1;
  }
  nbrOfHeaderLines <- kk;
  close(fh);

  fh <- file(filename, "r");
  nbrOfColumns <- length(header);
  df <- scan(file=fh, what=rep(list(""), nbrOfColumns+1), sep="\t", quote="", quiet=TRUE, strip.white=TRUE, blank.lines.skip=TRUE, skip=nbrOfHeaderLines);
  close(fh);

  # Remove the first column which contains all "SPOT".
  df <- df[-1];

  # Convert all columns.
  for (k in 1:nbrOfColumns) {
    x <- type.convert(as.character(df[[k]]), na.strings=c("Inf","NA"));
    x[df[[k]]=="Inf"] <- Inf;
    df[[k]] <- x;
  }

  # Add a 'slide' column
#  slide <- rep(1, length(df[[1]]));
#  df <- c(list(slide), df);

  # Put names on columns
  names(df) <- c(header);
#  names(df) <- c("slide", header);
  
  list(header=header, remark=remark, spot=as.data.frame(df), version="1.0");
}, private=TRUE, static=TRUE);



#########################################################################/**
# @RdocMethod read
#
# @title "Reads one or several ScanAlyze files into a ScanAlyzeData object"
#
# \description{
#  @get "title".
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
# \value{Returns a @see "ScanAlyzeData" object.}
#
# @author
#
# \examples{
#   sa <- ScanAlyzeData$read(pattern="group.*.dat", path=system.file("data-ex", package="aroma"))
#   raw <- getRawData(sa)
# }
#
# \references{
#   The two data files in the example above was contributed by
#   Michael Stadler, Bioinformatics Group, Swiss Institute for
#   Experimental Cancer Research.
# }
#
# \seealso{
#   For pattern formats see @see "base::list.files".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("read", "ScanAlyzeData", function(static, filename=NULL, path=NULL, pattern=NULL, verbose=FALSE, ...) {
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
  
  # Support both path as a single string or as a vector of character strings.
  path <- rep(path, length.out=length(filename));

  res <- NULL;
  for (k in seq(length(filename))) {
    gc(); # Call the garbage collector in case we are running low in memory.
    tmp <- readOneFile(static, filename=filename[k], path=path[k], verbose=verbose, ...);
    slidename <- basename(filename[k]);
    setSlideNames(tmp, slidename);
    if (is.null(res))
      res <- tmp
    else {
      append(res, tmp);
    }
  }

  gc();

  res;
}, static=TRUE);


# Kept for backward compatibility.
setMethodS3("readAll", "ScanAlyzeData", function(...) {
  ScanAlyzeData$read(...);
}, private=TRUE, static=TRUE, deprecated=TRUE)




#########################################################################/**
# @RdocMethod getRawData
#
# @title "Gets the raw (foreground and background) intensites"
#
# \description{
#  Extracts the red and green spot intensitites (both foreground and background)
#  and returns a @see "RawData" object.
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{Specifying which slides to be extracted. If @NULL, 
#     all slides are considered.}
#   \item{fg}{If \code{"auto"}, the default foreground estimates according
#     to @seemethod "getForeground" is returned.}
#   \item{bg}{If \code{"auto"}, the default foreground estimates according
#     to @seemethod "getBackground" is returned.}
# }
#
# \value{
#   Returns a @see "RawData" object containing the specified slides.
# }
#
# \details{
#   The R and Rb channels will come from the CH2* fields, and
#   the G and Gb channels will come from the CH1* fields.
#   To swap the channels just use dyeSwap().
# }
#
# \examples{
#   sa <- ScanAlyzeData$read("group4.dat", path=system.file("data-ex", package="aroma"))
#
#   # Get the foreground and the background
#   raw <- getRawData(sa)
#
#   # The the background corrected data
#   ma <- getSignal(raw, bgSubtract=FALSE)
#
#   # Plot M vs A with a lowess line through the data points
#   plot(ma, slide=1)
#   lowessCurve(ma, lwd=2, gridwise=TRUE)
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getRawData", "ScanAlyzeData", function(this, slides=NULL, fg=c("auto"), bg=("auto")) {
  rgFg <- getForeground(this, which=fg, slides=slides);
  rgBg <- getBackground(this, which=bg, slides=slides);
  
  RawData(R=rgFg$R, G=rgFg$G, Rb=rgBg$R, Gb=rgBg$G, layout=getLayout(this), 
                                                       extras=this$.extras);
});


setMethodS3("as.RawData", "ScanAlyzeData", function(this, ...) {
  getRawData(this, ...);
})


#########################################################################/**
# @RdocMethod plotSpatial
#
# @title "Creates a spatial plot of a slide"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{slide}{The slide to be plotted.}
#   \item{pch}{The spot symbol. Default is \code{20} (solid disk).}
#   \item{col}{The color of the spots. If @NULL, default color is used.}
#   \item{palette}{If \code{col} is not specified, colors are generated
#   automaticially from the signals in the two foreground channels.
#   If \code{redgreen}, a red to green colors scale will be used.
#   If \code{blueyellow}, a blue to yellow colors scale will be used.}
#   \item{A.range}{The range of the log (base 2) spot intensities. Used
#    only for generating colors.}
#   \item{M.range}{The range of the log (base 2) spot ratios. Used
#    only for generating colors.}
#   \item{yaxt, ...}{Parameters accepted by \code{plot}.}
# }
#
# \value{Returns nothing.}
#
# \details{
#   The ScanAlyze software does not return the center position of a spot,
#   but the \code{TOP} \code{LEFT} \code{BOT} and \code{RIGHT} coordinates.
#   This method assumes that the center of the spot is in the center of
#   this box.
# }
#
# \examples{
#   sa <- ScanAlyzeData$read("group4.dat", path=system.file("data-ex", package="aroma"))
#
#   subplots(2)
#
#   opar <- par(bg="black")
#   plotSpatial(sa)
#   par(opar)
#
#   opar <- par(bg="black")
#   plotSpatial(sa, palette="blueyellow")
#   par(opar)
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("plotSpatial", "ScanAlyzeData", function(this, slide=1, pch=20, yaxt=NULL, col=NULL, palette="redgreen", A.range=c(0,16), M.range=c(-1,1), ...) {
  slide <- validateArgumentSlides(this, slide=slide);

  if (missing(yaxt))
    yaxt0 <- "n";
  
  if (is.null(col)) {
    rg <- getForeground(this, slides=slide);
    R <- rg$R[,1];
    G <- rg$G[,1];
    rm(rg);
    col <- MicroarrayData$createColors(R,G, type="RG", palette=palette, A.range=A.range, M.range=M.range);
  }

  xy <- getSpotPosition(this, slides=slide);  
  # Upper *left* corner is (0,0)
  x <- xy$x; y <- -xy$y; # Will also be used as xlab and ylab.

  plot(x,y, pch=pch, yaxt=yaxt0, col=col, ...);
  
  if (missing(yaxt)) {
    # Add the y-axis with the correct ticks
    yaxp <- par("yaxp");
    yaxis <- seq(yaxp[1],yaxp[2], length=yaxp[3]+1);
    axis(2, at=yaxis, labels=-yaxis);
  }
})


setMethodS3("plotSpatial3d", "ScanAlyzeData", function(this, field="RAT2", ...) {
  NextMethod("plotSpatial3d", this, field=field, ...);
})



setMethodS3("normalizeGenewise", "ScanAlyzeData", function(this, fields=NULL, bias=0, scale=1, ...) {
  if (is.null(fields)) {
    fields <- getFieldNames(this);
    exclude <- c("SPOT", "GRID", "TOP", "LEFT", "BOT", "RIGHT", "ROW", "COL", "FLAG");
    fields <- setdiff(fields, exclude);
  }
  NextMethod("normalizeGenewise", this, fields=fields, bias=bias, scale=scale, ...);
})



############################################################################
#
#   S I G N A L    E S T I M A T E S
#
############################################################################
setMethodS3("getForeground", "ScanAlyzeData", function(this, which=c("mean","auto"), slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  which <- match.arg(which);
  if (which == "auto")
    which <- "mean";

  if (which == "mean") {
    fields <- c("CH2I", "CH1I");
  }

  # Assert that the fields do really exist.
  if (!all(fields %in% getFields(this))) {
    throw("The foreground estimates ", paste(fields, collapse=" and "), 
          " is not part of this ", data.class(this), " object.");
  }

  RGData(R=this[[fields[1]]][,slides], G=this[[fields[2]]][,slides], layout=getLayout(this));
})


# Default value is \code{"median"} as recommended to be used by the 
# author of the ScanAlyze software.
setMethodS3("getBackground", "ScanAlyzeData", function(this, which=c("median", "mean", "auto"), slides=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  which <- match.arg(which);
  if (which == "auto")
    which <- "median";

  if (which == "mean") {
    fields <- c("CH2BA", "CH1BA");
  } else if (which == "median") {
    fields <- c("CH2B", "CH1B");
  }

  # Assert that the fields do really exist.
  if (!all(fields %in% getFields(this))) {
    throw("The background estimates ", paste(fields, collapse=" and "), 
          " is not part of this ", data.class(this), " object.");
  }

  RGData(R=this[[fields[1]]][,slides], G=this[[fields[2]]][,slides], layout=getLayout(this));
})



############################################################################
#
#   S P O T   G E O M E T R Y
#
############################################################################
setMethodS3("getArea", "ScanAlyzeData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  # "'SPIX' - Number of pixels in the foreground, i.e. in the spot."
  as.matrix(this[["SPIX"]][include,slides]);
})



setMethodS3("getBgArea", "ScanAlyzeData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  # "'BGPIX' - Number of pixels in the background."
  as.matrix(this[["BGPIX"]][include,slides]);
})


setMethodS3("getDiameter", "ScanAlyzeData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  # "'TOP, LEFT, BOT, RIGHT' - Coordinates of the box containing spot
  #                            ellipse, in image coordinates."
  # (Warning! This is a very bad approximation of the diameter of the
  #  spots, especially when fixed size boxed are used.)
  top    <- this[["TOP"]][include,slides];
  left   <- this[["LEFT"]][include,slides];
  bottom <- this[["BOT"]][include,slides];
  right  <- this[["RIGHT"]][include,slides];
  diameter <- ((right-left) + (top-bottom)) / 2;

  diameter;
})

setMethodS3("getCircularity", "ScanAlyzeData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  area <- getArea(this, slides=slides)[include];
  d <- getDiameter(this, slides=slides)[include];
  
  # area = pi*r^2 = pi*(d/2)^2 = pi/4 * d^2  <=> pi/4 = area / d^2
  # However, this is only true for circles. So calculate
  #   circularity = area / (pi/4 * d^2)
  # where 0 <= circularity <= 1, but for discretization reasons it
  # might be larger than one too.
  circularity <- area / (pi/4 * d^2)
  as.matrix(circularity)
})




setMethodS3("getGrid", "ScanAlyzeData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  this[["GRID"]][include,slides];
})

setMethodS3("getSpotRow", "ScanAlyzeData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  this[["ROW"]][include,slides];
})

setMethodS3("getSpotColumn", "ScanAlyzeData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  this[["COL"]][include,slides];
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
# \details{
#   The ScanAlyze software does not return the center position of a spot,
#   but the \code{TOP} \code{LEFT} \code{BOT} and \code{RIGHT} coordinates.
#   This method assumes that the center of the spot is in the center of
#   this box.
# }
#
# \examples{
#   sa <- ScanAlyzeData$read("group4.dat", path=system.file("data-ex", package="aroma"))
#
#   # Gets the spot positions
#   xy <- getSpotPosition(sa)
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getSpotPosition", "ScanAlyzeData", function(this, slides=NULL, index=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (is.null(index)) {
    index <- 1:nbrOfSpots(this);
  } else if (any(index < 1) || any(index > nbrOfSpots(this))) {
    throw("Argument 'index' is out of range.");
  }
  
  x <- (this[["LEFT"]][index,slides]+this[["RIGHT"]][index,slides]) / 2;
  y <- (this[["TOP"]][index,slides]+this[["BOT"]][index,slides]) / 2;
  
  SpotPosition(x=x, y=y)
})




############################################################################
# HISTORY:
# 2008-01-15
# o Replaced obsolete trycatch() with tryCatch().
# o BUG FIX: getDiameter() of ScanAlyzeData used a non-existing variable.
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
# 2003-10-27
# o Updated the algorithm to identify the number of grid rows. The former
#   used the position of the *last* spot in each grid, but in ScanAlyze
#   it is allowed to have missing spots and if the last printdips are non-
#   existing, the position will be NAs. Better to use the position of the
#   *first* spot instead.
# o Made the code more modularized so it can be reused by future versions
#   such as ScanAlyze20Data.
# o Added ScanAlyze20Data to be able to read ScanAlyze v2.0 data too.
#   It extends this class.
# 2003-10-26
# o Added support to read *.gz files.
# 2003-10-19
# o getBackground() now returns the *median* background estimates.
# 2003-06-15
# o TYPO FIX: Updated Rdoc about 'CH*KSD' and 'CH*KSP'.
# 2002-12-21
# o Added plotSpatial3d().
# o When reading slides from files each slide is now named as the filename.
# 2002-12-05
# o BUG FIX: The get*() methods did not return a matrix if there was only
#   one slide.
# 2002-09-24
# o Changed the attribute 'path="."' to 'path=""' in read().
# o Update the Rdoc comments for as.RawData() to be more explicit about
#   which channel is used for R and which is used for G.
# 2002-09-21
# o read() can now read one or several files specified by names or by a
#   pattern. This is identical to readAll(), which is now just calling
#   read() for backward compatibilities.
# 2002-08-20
# o Replace 'append(super(this),other)' with NextMethod("append") in
#   method append().
# 2002-06-24
# o BUG FIX: getSpotPosition() was broken.
# 2002-05-28
# o Added some Rdoc comments on the CH1EDGEA, CH2EDGEA fields.
# 2002-05-03
# o Added getArea(), getCircularity() and getBgArea().
# 2002-04-21
# * Added question marks to the unknown fields in the Rdoc for this class.
# * Added trial version of normalizeGenewise().
# 2002-04-03
# * read() and write() now supports File objects.
# 2002-03-25
# * Added getSpotPositions() and plotSpatial().
# 2002-02-28
# * Made class implements interface Morpable.
# 2002-02-26
# * Added example data for ScanAlyze too.
# * Added write().
# * Updated the code and the Rdoc's to make use of setMethodS3().
# 2001-11-17
# * Updated readAll() to also include pattern matching.
# 2001-11-13
# * Rdoc bug: Wrote CH2B and CH2BA instead of CH1B and CH1BA.
# 2001-11-12
# * Added readAll.ScanAlyzeData().
# 2001-07-12
# * Bug fix: Forgot to do quote="" in all scan() and read.table() calls.
#   Trying to read files with cells including "'" would give an error.
# 2001-07-10
# * Create a safe and a fast read method. The safe is too slow!
# * Created from ScanAlyzeData. Javier needed it.
#
############################################################################
