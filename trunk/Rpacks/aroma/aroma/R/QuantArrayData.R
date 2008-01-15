#########################################################################/**
# @RdocClass QuantArrayData
#
# @title "The QuantArrayData class"
#
# @synopsis
#
# \description{
#  @classhierarchy
#
#  Creates an empty QuantArrayData object.
# }
#
# \arguments{
#   \item{layout}{A @see "Layout" object specifying the spot layout of the
#    slides in this data set.}
# }
#
# \section{Fields}{
#   \item{ch1 Intensity, ch2 Intensity, ...}{The spot intensity.
#     What the authors mean by intensity is not clear, but on p46 in [2]
#     it says "Quantitation. Output allows the user to choose the way
#     the spot intensities are calculated: a) the sum of the intensities
#     for the various pixels of the spot and the background (Total
#     Intensities), b) the average of the intensities (Mean Intensity),
#     c) the most frequent occurrence of a pixel value (Mode Intensity),
#     or d) the Median Intensity can be reported."
#
#     If this is true, then where is the information about which method
#     was used? Looking at the header of a typical QuantArray data file
#     this information is not stored anywhere.
#   }
#   \item{ch1 Intensity Std Dev, ch2 Intensity Std Dev, ...}{
#     The standard deviation of the spot pixel intensities.}
#   \item{ch1 Spot Uniformity, ch2 Spot Uniformity, ...}{
#     Value in [0,1].
#     The 0.99 quantile minus the 0.01 quantile pixel signal divided
#     by the pixel intensity range.}
#
#   \item{ch1 Background, ch2 Background, ...}{The background intensity.}
#   \item{ch1 Background Std Dev, ch2 Background Std Dev, ...}{
#     The standard deviation of the background pixel intensities.}
#   \item{ch1 Bkg. Uniformity, ch2 Bkg. Uniformity, ...}{Value in [0,1].}
#
#   \item{ch1 Diameter, ch2 Diameter, ...}{The average spot diameter
#     (in micrometers).}
#   \item{ch1 Area, ch2 Area, ...}{The spot area as the number of pixels
#     (in square micrometers).}
#   \item{ch1 Footprint, ch2 Footprint, ...}{The distance (in micrometers)
#     between the expected position of the spot and its actual position.}
#   \item{ch1 Circularity, ch2 Circularity, ...}{
#     Value in [0,1]. The circularity.}
#   \item{ch1 Replicate Uniformity, ch2 Replicate Uniformity, ...}{
#     Value in [0,1]. Uniformity of replicated spots.}
#
#   \item{ch1 Confident, ch2 Confident, ...}{
#     Value in [0,1].
#     Different methods may have been used. See in
#     \code{sections[["Protocol Info"]]} of the QuantArrayData object
#     which was used.
#     Let \code{(c[1], c[2], ..., c[Q])} be the individual measurement
#     confidences. Then
#     "Minimum") \code{cAll = min(c)}.
#     "Weighted Average") \code{cAll = sum(w*c)/sum(w)} where
#     \code{w} are weights.
#     "Product") \code{cAll = (prod(c))^(1/Q)}.
#   }
#
#  \emph{Redundant fields}:\cr
#   \item{ch1 Percentage, ch2 Percentage, ...}{\bold{REDUNDANT!}
#     \code{=="ch1 Intensity"/("ch1 Intensity"+"ch2 Intensity")} and so on.
#     For each spot, the percentage of the total brightness (the sum of all
#     channels' intensities) for each channel.}
#   \item{ch1 Ratio, ch2 Ratio, ...}{\bold{REDUNDANT!}
#     \code{=="ch2 Intensity"/"ch1 Intensity"} or similar.
#     For each spot, the relative intensity of a spot, compared to the
#     corresponding spot in the control. The ratio value for the channel
#     that is the control channel (commonly channel 1) are all equal to 1.
#     The ratio value for the other channels are equal to their relative
#     ratio of the intensity compared to the control.
#     Background subtraction might have been applied first, but that is
#     still reduandant.}
#   \item{ch1 Signal Noise Ratio, ch2 Signal Noise Ratio, ...}{
#     \bold{REDUNDANT!}
#     \code{=="ch1 Intensity"/"ch1 Background Std Dev"} and so on.
#     The ratio of the spot intensity to the standard deviation of the
#     local background of all spots in the microarray.}

#   Sources: Our own "research" and [2].
# }

# \section{Methods}{
#  @allmethods
# }
#
# @author
#
# \references{
#   [1] QuantArray software by PerkinElmer Life Sciences, 
#       \url{http://lifesciences.perkinelmer.com/}.\cr
#   [2] Packard BioScience,
#       QuantArray Microarray Analysis Software Manual, 2001.\cr
# }
#
# \examples{
#  \dontrun{
#   # At the moment there is no QuantArray sample files in the package...
#   # qa <- QuantArrayData$read("quantarray123.txt", path=system.file("data-ex", package="aroma"))
#   # ...will use a GenePix sample file instead.
#   qa <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#
#   # Get the raw data; Rfg, Gfg, Rbg and Gbg
#   raw <- getRawData(qa)
#
#   # The the background corrected data
#   ma <- getSignal(raw)
#
#   # Plot M vs A with a lowess line through the data points
#   subplots(4)
#   plot(ma)
#   plotSpatial(ma)
#   boxplot(ma, groupBy="printtip")
#  }
# }
#
# \note{
#   The laser beam with wavelength 635nm is the red laser and excites the
#   Cy5 dye, and the one with wavelength 532nm is the green laser which
#   excites the Cy3 dye.
# }
#*/#########################################################################
setConstructorS3("QuantArrayData", function(layout=NULL) {
  extend(MicroarrayData(layout=layout), "QuantArrayData",
    version  = list(),
    header   = list(),
    sections = list()
  )
})




setMethodS3("read.table", "QuantArrayData", function(this, file, header=FALSE, sep="", quote="\"'", dec=".", row.names, col.names, as.is=FALSE, na.strings="NA", colClasses=NA, nrows=-1, skip=0, check.names=TRUE, fill=!blank.lines.skip, strip.white=FALSE, blank.lines.skip=TRUE, comment.char="#", flush=FALSE) {
  # Open the file
  if (is.character(file)) {
    file <- file(file, "r")
    on.exit(close(file))
  }
  
  if (!inherits(file, "connection")) 
    throw("Argument `file' must be a character string or connection.")
  if (!isOpen(file)) {
    open(file, "r")
    on.exit(close(file))
  }

  # Skip initial rows
  if (skip > 0) 
    readLines(file, skip)

  # Get the number of columns from the first 5 rows...
  if (R.Version()$major <= 1 && R.Version()$minor < 9.0) {
    lines <- .Internal(readTableHead(file, 5, comment.char, blank.lines.skip));
  } else if (R.Version()$major == 2 && R.Version()$minor < 1.0) {
    lines <- .Internal(readTableHead(file, 5, comment.char, blank.lines.skip));
  } else {
    lines <- .Internal(readTableHead(file, 5, comment.char, blank.lines.skip, quote, sep));
  }
  nlines <- length(lines)
  if (!nlines) {
    if (missing(col.names)) 
      throw("No lines available in input.")
    else {
      tmp <- vector("list", length(col.names))
      names(tmp) <- col.names
      class(tmp) <- "data.frame"
      return(tmp)
    }
  }
  if (all(nchar(lines) == 0)) 
    throw("Empty beginning of file.")
  pushBack(c(lines, lines), file)
  first <- scan(file, what = "", sep = sep, quote = quote, nlines = 1, quiet = TRUE, skip = 0, strip.white = TRUE, blank.lines.skip = blank.lines.skip, comment.char = comment.char)
  col1 <- if (missing(col.names)) 
    length(first)
  else
    length(col.names)
  col <- numeric(nlines - 1)
  for (i in seq(along = col))
    col[i] <- length(scan(file, what = "", sep = sep, quote = quote, nlines = 1, quiet = TRUE, skip = 0, strip.white = strip.white, blank.lines.skip = blank.lines.skip, comment.char = comment.char))
  cols <- max(col1, col)

  rlabp <- (cols - col1) == 1
  if (rlabp && missing(header)) 
    header <- TRUE
  if (!header) 
    rlabp <- FALSE
  if (header) {
    readLines(file, 1)
    if (missing(col.names)) 
      col.names <- first
    else if (length(first) != length(col.names)) 
      warning("header and `col.names' are of different lengths")
  }
  else if (missing(col.names)) 
    col.names <- paste("V", 1:cols, sep = "")
  if (length(col.names) + rlabp < cols) 
    throw("More columns than column names.")
  if (fill && length(col.names) > cols) 
    cols <- length(col.names)
  if (!fill && cols > 0 && length(col.names) > cols) 
    throw("More column names than columns.")
  if (cols == 0) 
    throw("First five rows are empty: giving up.")
  if (check.names) 
    col.names <- make.names(col.names)
  if (rlabp) 
    col.names <- c("row.names", col.names)
  if (length(colClasses) < cols) 
    colClasses <- rep(colClasses, len = cols)
  what <- rep(list(""), cols)
  names(what) <- col.names
  colClasses[colClasses %in% c("real", "double")] <- "numeric"
  known <- colClasses %in% c("logical", "integer", "numeric", "complex", "character")
  what[known] <- sapply(colClasses[known], do.call, list(0))

  data <- scan(file = file, what = what, sep = sep, quote = quote, dec = dec, nmax = nrows, skip = 0, na.strings = na.strings, quiet = TRUE, fill = fill, strip.white = strip.white, blank.lines.skip = blank.lines.skip, multi.line = FALSE, comment.char = comment.char, flush=flush)
  
  nlines <- length(data[[1]])
  if (cols != length(data)) {
    warning(paste("cols =", cols, " != length(data) =", length(data)))
    cols <- length(data)
  }
  if (is.logical(as.is)) {
    as.is <- rep(as.is, length = cols)
  }
  else if (is.numeric(as.is)) {
    if (any(as.is < 1 | as.is > cols)) 
      throw("Invalid numeric as.is expression.")
    i <- rep(FALSE, cols)
    i[as.is] <- TRUE
    as.is <- i
  }
  else if (length(as.is) != cols) 
    throw("as.is has the wrong length: ", length(as.is), "!= cols =", cols);
  for (i in 1:cols) {
    if (known[i]) 
      next
    data[[i]] <- if (!is.na(colClasses[i])) 
      as(data[[i]], colClasses[i])
    else type.convert(data[[i]], as.is = as.is[i], dec = dec)
  }
  if (missing(row.names)) {
    if (rlabp) {
      row.names <- data[[1]]
      data <- data[-1]
    } else
      row.names <- as.character(seq(len = nlines))
  }
  else if (is.null(row.names)) {
    row.names <- as.character(seq(len = nlines))
  }
  else if (is.character(row.names)) {
    if (length(row.names) == 1) {
       rowvar <- (1:cols)[match(col.names, row.names, 0) == 1]
       row.names <- data[[rowvar]]
       data <- data[-rowvar]
    }
  }
  else if (is.numeric(row.names) && length(row.names) == 1) {
    rlabp <- row.names
    row.names <- data[[rlabp]]
    data <- data[-rlabp]
  }
  else throw("Invalid row.names specification.")
  class(data) <- "data.frame"
  row.names(data) <- row.names
  data
}, protected=TRUE, static=TRUE) # read.table()



                                                     
#########################################################################/**
# @RdocMethod readOneFile
#
# @title "Reads a QuantArray result file"
#
# \description{
#  @get "title".
#  Currently QuantArray v2 and v3 are supported. If the version is not
#  recognized it will try to make a best guess and read it in anyway.
#  However, it can not be guaranteed that all fields will be of the correct
#  data type.
#
#  This method also reads Unicoded QuantArray files. If the file is
#  Unicoded it is first translated into a temporary ASCII file, which is
#  then read. The translation from Unicode to ASCII is done by "brute force",
#  i.e. by excluding the highbyte and only keeping the lowbytes. This means
#  that some characters will be incorrectly translated. If that happens, a
#  warning will be given, otherwise not.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file to be read.}
#   \item{path}{The path to the file.}
#   \item{verbose}{If @TRUE, information will printed out during
#                  the reading of the file.}
# }
#
# \value{Returns a @see "QuantArrayData" object.}
#
# @author
#
# \examples{
#  \dontrun{
#   # At the moment there is no QuantArray sample files in the package...
#   # qa <- QuantArrayData$read("quantarray123.txt", path=system.file("data-ex", package="aroma"))
#   # ...will use a GenePix sample file instead.
#   qa <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#
#   raw <- getRawData(qa)
#  }
# }
#

# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("readOneFile", "QuantArrayData", function(this, filename, path=NULL, verbose=FALSE) {
  # Unicode references:
  #  [1] http://www.cl.cam.ac.uk/~mgk25/unicode.html
  isUnicoded <- function(filename) {
    bfr <- readBin(filename, what=integer(), size=2, n=1, endian="little");
    res <- TRUE;
    if (bfr == -257)     # 0xFFFE
      attr(res, "endian") <- "little"
    else if (bfr == -2)  # 0xFEFF
      attr(res, "endian") <- "big"
    else
      res <- FALSE;
    return(res);
  }
  
  readHeader <- function(filename, beginRow, endRow, sep="\t", quote="") {
    skip <- beginRow;
    nlines <- endRow-beginRow-1;
    raw <- scan(filename, what=character(0), skip=skip, nlines=nlines, sep="\n", quote=quote, quiet=TRUE);
    raw <- as.character(raw);
    raw <- strsplit(raw, split=sep);
    names <- unlist(lapply(raw, FUN=function(x) x[1]));
    pinfo <- lapply(raw, FUN=function(x) x[-1]);
    names(pinfo) <- names;
    pinfo;
  }

  
  readSection <- function(filename, beginRow, endRow, sep="\t", quote="") {
    skip <- beginRow;
    nlines <- endRow-beginRow-1;
    raw <- scan(filename, what=character(0), skip=skip, nlines=nlines, sep="\n", quote=quote, quiet=TRUE);
    raw <- as.character(raw);
    raw <- strsplit(raw, split=sep);
    raw;
  }
  

  # The 'Measurements' section is tricky since some of the last rows normally contains
  # several trailing TAB's. Why this is I do not know, but it makes read.table() to
  # get confused and complain.
  readMeasurements <- function(filename, beginRow, endRow, sep="\t", quote="") {
    field <- list("Number"="integer", "Array Row"="integer", "Array Column"="integer", "Row"="integer", "Column"="integer", "Name"="character", "ch1 Ratio"="double", "ch1 Percent"="double", "ch2 Ratio"="double", "ch2 Percent"="double", "Ignore Filter"="integer", "NA"=NA);

    # Read the header line
    skip <- beginRow;
    df <- read.table(filename, header=TRUE, skip=skip, nrows=1, sep=sep, quote=quote, check.names=FALSE);
    header <- names(df);
  
    unknown <- which(is.na(match(header, names(field))));
    if (length(unknown) > 0) {
      warning(paste("Unknown 'Measurement' field(s) found in QuantArray file. Will try to make an intelligent guess about the data type: ", paste(header[unknown], collapse=", "), sep=""));
      header[unknown] <- "NA";
    }

#    header <- c(header, rep("NA", 10));
    colClasses <- unlist(field[header]);
    nrows <- endRow-skip-2;
    df <- QuantArrayData$read.table(filename, colClasses=colClasses, header=TRUE, skip=skip, nrows=nrows, sep=sep, quote=quote, check.names=FALSE, fill=FALSE, flush=TRUE);  # fill=FALSE will force errors.
#    colnames(df) <- gsub("[.]", " ", colnames(df));
    df;
  } # readMeasurements()


  readData <- function(filename, beginRow, endRow, sep="\t", quote="") {
    # QuantArray v2 fields
#    field <- list("Number"="integer", "Array Row"="integer", "Array Column"="integer", "Row"="integer", "Column"="integer", "Name"="character", "X Location"="integer", "Y Location"="integer", "ch1 Intensity"="double", "ch1 Background"="double", "ch1 Intensity Std Dev"="double", "ch1 Background Std Dev"="double", "ch2 Intensity"="double", "ch2 Background"="double", "ch2 Intensity Std Dev"="double", "ch2 Background Std Dev"="double", "Ignore Filter"="integer", "NA"=NA);
    field <- list("Number"="integer", "Array Row"="integer", "Array Column"="integer", "Row"="integer", "Column"="integer", "Name"="character", "X Location"="double", "Y Location"="double", "ch1 Intensity"="double", "ch1 Background"="double", "ch1 Intensity Std Dev"="double", "ch1 Background Std Dev"="double", "ch2 Intensity"="double", "ch2 Background"="double", "ch2 Intensity Std Dev"="double", "ch2 Background Std Dev"="double", "Ignore Filter"="integer", "NA"=NA);

    # QuantArray v3 fields
#    field <- c(field, "ch1 Diameter"="double", "ch1 Area"="integer", "ch1 Footprint"="double", "ch1 Circularity"="double", "ch1 Spot Uniformity"="double", "ch1 Bkg. Uniformity"="double", "ch1 Signal Noise Ratio"="double", "ch1 Confidence"="double", "ch2 Diameter"="double", "ch2 Area"="integer", "ch2 Footprint"="double", "ch2 Circularity"="double", "ch2 Spot Uniformity"="double", "ch2 Bkg. Uniformity"="double", "ch2 Signal Noise Ratio"="double", "ch2 Confidence"="double");
    field <- c(field, "ch1 Diameter"="double", "ch1 Area"="double", "ch1 Footprint"="double", "ch1 Circularity"="double", "ch1 Spot Uniformity"="double", "ch1 Bkg. Uniformity"="double", "ch1 Signal Noise Ratio"="double", "ch1 Confidence"="double", "ch2 Diameter"="double", "ch2 Area"="double", "ch2 Footprint"="double", "ch2 Circularity"="double", "ch2 Spot Uniformity"="double", "ch2 Bkg. Uniformity"="double", "ch2 Signal Noise Ratio"="double", "ch2 Confidence"="double");

    # Read the header line
    skip <- beginRow;
    df <- read.table(filename, header=TRUE, skip=skip, nrows=1, sep=sep, quote=quote, check.names=FALSE);
    header <- names(df);
  
    unknown <- which(is.na(match(header, names(field))));
    if (length(unknown) > 0) {
      warning(paste("Unknown 'Data' field(s) found in QuantArray file. Will try to make an intelligent guess about the data type: ", paste(header[unknown], collapse=", "), sep=""));
      header[unknown] <- "NA";
    }

    colClasses <- unlist(field[header]);
    nrows <- endRow-skip-2;
    df <- QuantArrayData$read.table(filename, colClasses=colClasses, header=TRUE, skip=skip, nrows=nrows, sep=sep, quote=quote, check.names=FALSE, flush=TRUE);
#    colnames(df) <- gsub("[.]", " ", colnames(df));
    df;
  } # readData()

  filename <- Arguments$getReadablePathname(filename, path);  

  # Support gzip'ed files too.
  if (regexpr("[.]gz$", filename) != -1) {
    tmpname <- tempfile();
    n <- gunzip(filename, tmpname);
    filename <- tmpname;
    on.exit(file.remove(tmpname));
  }
  
  if (verbose) cat("Reading file ", filename, "...\n", sep="");

  isUnicoded <- isUnicoded(filename);
  if (isUnicoded) {
    if (verbose) cat("Converting unicoded file to a temporary ASCII file...\n", sep="");
    len <- file.info(filename)$size;
    len <- len %/% 2;
    
    filework <- tempfile();
    on.exit(unlink(filework));
    fhIn <- file(filename, "rb");
    seek(fhIn, where=2);
    fhOut <- file(filework, "wb");
    chunkSize <- 65536;   # Too save memory!
    missTranslation <- FALSE;
    while (len > 0) {
      bfr <- readBin(fhIn, what=integer(), size=2, n=chunkSize, endian=attr(isUnicoded, "endian"));
      if (!missTranslation) {
        if (any(bfr %/% 256 != 0)) {
          missTranlation <- TRUE;
          warning(paste("The file", filename, "is unicoded, but [R] only support ASCII and some of the characters found in the unicoded file was by \"brute force\" (ignoring the highbyte) misstranslated into ASCII."));
        }
      }
      bfr <- as.integer(bfr %% 256);
      writeBin(bfr, con=fhOut, size=1);
      len <- len - chunkSize;
    }
    close(fhIn);
    close(fhOut);
    if (verbose) cat("Converting unicoded file to a temporary ASCII file...ok\n", sep="");
  } else {
    filework <- filename;
  }

  fh <- file(filework, "r");

  ########################################################################
  # Step 1 - Prescan the file for Begin and End tags
  ########################################################################
  # Read the whole file into a number of strings containing the rows.
  rows <- readLines(fh);
  close(fh);

  # Find all Begin tags
  pos <- regexpr("^[ \t]*Begin[ ]*", rows);
  beginRow <- which(pos != -1);
  rows[beginRow] <- gsub("^[ \t]+", "", rows[beginRow]);
  rows[beginRow] <- gsub("[ \t]+$", "", rows[beginRow]);
  attr(beginRow, "match.length") <- attr(pos, "match.length")[beginRow];
  
  # Find all End tags
  pos <- regexpr("^[ \t]*End[ ]*", rows);
  endRow <- which(pos != -1);
  rows[endRow] <- gsub("^[ \t]+", "", rows[endRow]);
  rows[endRow] <- gsub("[ \t]+$", "", rows[endRow]);
  attr(endRow, "match.length") <- attr(pos, "match.length")[endRow];
  
  # Verify that the Begin tags match the End tags.
  len <- length(beginRow);
  if (len != length(endRow))
    throw("QuantArray file format error: The number of Begin tags (", len, ") does not match the number of End tags (", length(endRow), ").");

  beginTag <- substring(rows[beginRow], attr(beginRow, "match.length")+1);
  endTag <- substring(rows[endRow], attr(endRow, "match.length")+1);

  if (length(union(beginTag, endTag)) != len)
    throw("QuantArray file format error: The Begin tags do not match the End tags.");
  
  # Don't need the 'rows' anymore.
  rm(rows); # Save memory.

  # Verify that the exists a "Data" Begin/End tag.
  if (!any(beginTag == "Data"))
    throw("QuantArray file format error: Expected a 'Data' section.");
  

  ########################################################################
  # Step 1b - Read header info
  ########################################################################
  header <- readSection(filework, 0, beginRow[1]-1);
  keys <- lapply(header, FUN=function(x) x[1]);
  header <- lapply(header, FUN=function(x) x[-1]);
  names(header) <- keys;
  if (verbose)
    cat("Identified version ", header[["Version"]], ".\n", sep="");
  
  ########################################################################
  # Step 2 - Read the "Data"
  ########################################################################
  tag <- "Data";
  begin <- beginRow[beginTag == tag];
  end   <- endRow[endTag == tag];
  nbrOfRows <- end-begin-2;
  
  if (verbose == TRUE)
    cat("Reading ", nbrOfRows, " rows of 'Data'...", sep="");
  df <- readData(filework, begin, end);
  if (verbose == TRUE)
    cat("done\n")

  ngrid.r <- max(df[,"Array Row"]);
  ngrid.c <- max(df[,"Array Column"]);
  nspot.r <- max(df[,"Row"]);
  nspot.c <- max(df[,"Column"]);
  nbrOfSpots <- ngrid.r*ngrid.c*nspot.r*nspot.c;

  ########################################################################
  # Step 3 - Read (optional) "Measurements"
  ########################################################################
  tag <- "Measurements";
  begin <- beginRow[beginTag == tag];
  end   <- endRow[endTag == tag];
  nbrOfRows <- end-begin-2;
  if (nbrOfSpots != nbrOfRows) {
    warning("QuantArray file format error: According to the 'Protocol Info' tag the cDNA microarray contains ", nbrOfSpots, " spots, but there are only ", nbrOfRows, " rows of data in the 'Measurements' section of the file. Will ignore the header information (anyway).");
  }
  
  if (verbose == TRUE)
    cat("Reading ", nbrOfRows, " rows of 'Measurements'...", sep="");
  dfMeasurements <- readMeasurements(filework, begin, end);
  if (verbose == TRUE)
    cat("done\n")
  keep <- setdiff(names(dfMeasurements), names(df));

  df <- cbind(df, dfMeasurements[,keep]);
  rm(dfMeasurements);
  gc();

  ########################################################################
  # Step 4 - Read other (optional) sections
  ########################################################################
  tags <- c("Image Info", "Filter");
  tags <- setdiff(beginTag, c("Data", "Measurements"));
  sections <- list();
  for (tag in tags) {
    begin <- beginRow[beginTag == tag];
    end   <- endRow[endTag == tag];
    if (verbose == TRUE)
      cat("Reading ", end-begin-1, " rows of '", tag, "'...", sep="");
    sections <- c(sections, list(readSection(filework, begin, end)));
    if (verbose == TRUE)
      cat("done\n")
  } # for (tag in tags)
  names(sections) <- tags;

  ########################################################################
  # Step 5 - Creating the QuantArrayData object
  ########################################################################
  id <- as.character(df[,"Name"]);  # Make sure the id's are not factors!
  layout <- Layout(ngrid.r,ngrid.c, nspot.r,nspot.c, id=id)
  rm(ngrid.r, ngrid.c, nspot.r, nspot.c, id);
  
  qa <- QuantArrayData(layout=layout)
  qa$version <- header[["Version"]];
  qa$header <- header;
  qa$sections <- sections;
  rm(header, sections);

  qa$.fieldNames <- names(df);
  for (field in names(df)) {
    qa[[field]] <- as.matrix(df[[field]]);
    df[[field]] <- NULL; # Save memory
  }
  rm(df); gc(); # To minimize memory usage!

  if (verbose) cat("Reading file ", filename, "...ok\n", sep="");
  
  qa;
}, private=TRUE, static=TRUE) # readOneFile()




#########################################################################/**
# @RdocMethod read
#
# @title "Reads one or several QuantArray files into one QuantArrayData object"
#
# \description{
#  @get "title".
#  Currently QuantArray v2 and v3 are supported. If the version is not
#  recognized it will try to make a best guess and read it in anyway.
#  However, it can not be guaranteed that all fields will be of the correct
#  data type.
#
#  This method also reads Unicoded QuantArray files. If the file is
#  Unicoded it is first translated into a temporary ASCII file, which is
#  then read. The translation from Unicode to ASCII is done by "brute force",
#  i.e. by excluding the hi-bytes and only keeping the lo-bytes. This means
#  that some characters will be incorrectly translated. If that happens, a
#  warning will be given, otherwise not.
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
# \value{Returns a @see "QuantArrayData" object.}
#
# @author
#
# \examples{
#  \dontrun{
#   # Read two QuantArray data files
#   # filenames <- c("quantarray123.gpr", "quantarray123.gpr");
#   # qa <- QuantArrayData$read(filenames, path=system.file("data-ex", package="aroma"))
#
#   # Do not have any QuantArray sample files, will examplify using
#   # the GenePix sample files instead...
#   filenames <- c("gpr123.gpr", "gpr123.gpr");
#   qa <- GenePixData$read(filenames, path=system.file("data-ex", package="aroma"))
#
#   # Create a RawData object from this QuantArrayData objects.
#   raw <- getRawData(qa)
#  }
# }
#
# \seealso{
#   To write a slide to a QuantArray Results file
#   see @seemethod "write".
#   For pattern formats see @see "base::list.files".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("read", "QuantArrayData", function(this, filename=NULL, path=NULL, pattern=NULL, verbose=FALSE, ...) {
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
  
  # Support both path as a single string or as a vector of strings.
  path <- rep(path, length.out=length(filename));

  res <- NULL;
  for (k in seq(length(filename))) {
    gc(); # Call the garbage collector in case we are running low in memory.
    tmp <- QuantArrayData$readOneFile(filename=filename[k], path=path[k], verbose=verbose, ...);
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
setMethodS3("readAll", "QuantArrayData", function(...) {
  QuantArrayData$read(...);
}, private=TRUE, static=TRUE, deprecated=TRUE)




#########################################################################/**
# @RdocMethod write
#
# @title "Write a QuantArray result data file"
#
# \description{
#  Writes the QuantArrayData object to a file using QuantArray file format.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file to be written.}
#   \item{path}{The path to the file.}
#   \item{slide}{An @integer specifying which slide to be written to file.
#     Currently, only one slide at the time can be written.}
#   \item{...}{Arguments passed to \code{write.table}.}
# }
#
# \value{Returns nothing.}
#
# @author
#
# \examples{
#  \dontrun{
#   # At the moment there is no QuantArray sample files in the package...
#   # qa <- QuantArrayData$read("quantarray123.txt", path=system.file("data-ex", package="aroma"))
#   # ...will use a GenePix sample file instead.
#   qa <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#
#   # Writes the QuantArrayData object to a file named "temp.txt" in a format
#   # that is as close as possible to the original format.
#   write(qa, "temp.txt", overwrite=TRUE)
#  }
# }
#
# \seealso{
#   To read a QuantArray result data file see @seemethod "read".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("write", "QuantArrayData", function(this, filename, path=NULL, slide=1, overwrite=FALSE, ..., digits=9, verbose=FALSE) {
  slide <- validateArgumentSlide(this, slide=slide);
  filename <- Arguments$getWritablePathname(filename, path, mustNotExist=!overwrite);  

  fh <- file(filename, "w");
  on.exit(close(fh));

  # Write Header
  df <- t(as.data.frame(this$header));
  rownames(df) <- names(this$header);
  write.table(df, file=fh, append=TRUE, col.names=FALSE, sep="\t", quote=FALSE, row.names=TRUE);
  
  # Write Sections
  sectionNames <- names(this$sections);
  names <- setdiff(names(this$sections), c("Filter"));
  for (name in names) {
    cat(file=fh, "\n");
    cat(file=fh, "Begin ", name, "\n", sep="");
    section <- this$sections[[name]];
    for (k in seq(section))
      cat(file=fh, section[[k]], "\n", sep="\t")
    cat(file=fh, "End ", name, "\n", sep="");
  }

  
  dfLayout <- as.data.frame(getLayout(this));
  cols <- c("gridRow", "gridColumn", "spotRow", "spotColumn", "id");
  dfLayout <- dfLayout[,cols];
  dfLayout <- cbind(seq(nrow(dfLayout)), dfLayout);
  names(dfLayout) <- c("Number", "Array Row", "Array Column", "Row", "Column", "Name");
  
  knownFields <- list("Number"="integer", "Array Row"="integer", "Array Column"="integer", "Row"="integer", "Column"="integer", "Name"="character", "X Location"="integer", "Y Location"="integer", "ch1 Intensity"="double", "ch1 Background"="double", "ch1 Intensity Std Dev"="double", "ch1 Background Std Dev"="double", "ch2 Intensity"="double", "ch2 Background"="double", "ch2 Intensity Std Dev"="double", "ch2 Background Std Dev"="double", "Ignore Filter"="integer", "NA"=NA);
  knownFields <- c(knownFields, list("Number"="integer", "Array Row"="integer", "Array Column"="integer", "Row"="integer", "Column"="integer", "Name"="character", "ch1 Ratio"="double", "ch1 Percent"="double", "ch2 Ratio"="double", "ch2 Percent"="double", "Ignore Filter"="integer", "NA"=NA));

  oopt <- options(digits=digits);

  # Write Measurements
  fields <- c("ch1 Ratio", "ch1 Percent", "ch2 Ratio", "ch2 Percent", "Ignore Filter");
  df <- dfLayout;  
  for (field in fields) {
    data <- unlist(extract(this, field));
    type <- knownFields[[field]];
    if (is.null(type))
      type <- NA
    else if (type == "character")
      data <- as.character(data)
    else if (type == "integer")
      data <- as.integer(data)
    else if (type == "double")
      data <- as.double(data);
    df <- cbind(df, data);
    gc();
  }
  names <- c(names(dfLayout), fields);
  cat(file=fh, "\n");
  cat(file=fh, "Begin Measurements\n", sep="");
  cat(file=fh, names, sep="\t");
  cat(file=fh, "\n");
  write.table(df, file=fh, append=TRUE, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);
  cat(file=fh, "End Measurements\n", sep="");

  # Write Data
  fields <- setdiff(getFieldNames(this), fields);
  fields <- setdiff(fields, names(dfLayout));
  fields <- c(fields, "Ignore Filter");
  rm(df); gc();

  df <- dfLayout;
  for (field in fields) {
    data <- unlist(extract(this, field));
    type <- knownFields[[field]];
    if (is.null(type))
      type <- NA
    else if (type == "character")
      data <- as.character(data)
    else if (type == "integer")
      data <- as.integer(data)
    else if (type == "double")
      data <- as.double(data);
    df <- cbind(df, data);
    gc();
  }
  names <- c(names(dfLayout), fields);
  cat(file=fh, "\n");
  cat(file=fh, "Begin Data\n", sep="");
  cat(file=fh, names, sep="\t");
  cat(file=fh, "\n");
  write.table(df, file=fh, append=TRUE, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE);
  cat(file=fh, "End Data\n", sep="");

  rm(df); gc();

  options(oopt);
  
  # Write Sections
  sectionNames <- names(this$sections);
  names <- c("Filter");
  for (name in names) {
    cat(file=fh, "\n");
    cat(file=fh, "Begin ", name, "\n", sep="");
    section <- this$sections[[name]];
    for (k in seq(section))
      cat(file=fh, section[[k]], "\n", sep="\t")
    cat(file=fh, "End ", name, "\n", sep="");
  }
}) # write()





#########################################################################/**
# @RdocMethod getRawData
#
# @title "Gets the raw intensites from the QuantArray data structure"
#
# \description{
#  Extracts the red and green spot intensitites (both foreground and background)
#  from the QuantArrayData object and returns a @see "RawData" object.
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
#     If \code{"median"}, the median background intensities are returned.}
# }
#
# \value{
#   Returns a @see "RawData" object containing the specified slides.
# }
#
# \details{
#   The R and Rb channels will come from the ch1* fields, and
#   the G and Gb channels will come from the ch2* fields.
#   To swap the channels just use dyeSwap().
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getRawData", "QuantArrayData", function(qa, slides=NULL, fg="mean", bg="mean") {
  slides <- validateArgumentSlides(qa, slides=slides);
  if (!inherits(qa, "QuantArrayData"))
    throw("Argument 'qa' is not of class QuantArrayData: ", data.class(qa));

  #-------------------------------------------------------------------------
  # RGs as in the sma library.
  #-------------------------------------------------------------------------
  if (fg == "mean")
    fgs <- c("ch1 Intensity", "ch2 Intensity")
  else
    throw("Argument 'fg' has an unknown value: ", fg);

  if (bg == "mean")
    bgs <- c("ch1 Background", "ch2 Background")
  else
    throw("Argument 'bg' has an unknown value: ", bg);

  nas <- which(is.na(match(c(fgs,bgs), getFieldNames(qa))));
  if (length(nas) > 0) {
    throw("Strange, some fields do not exists in the QuantArrayData object 'qa'. Please, report this error to the author of the package: ", c(fgs,bgs)[nas]);
  }

  R  <- qa[[fgs[1]]][,slides];
  G  <- qa[[fgs[2]]][,slides];
  Rb <- qa[[bgs[1]]][,slides];
  Gb <- qa[[bgs[2]]][,slides];
    
  RawData(R=R, G=G, Rb=Rb, Gb=Gb, layout=getLayout(qa), extras=qa$.extras)
});


setMethodS3("as.RawData", "QuantArrayData", function(this, ...) {
  getRawData(this, ...);
})



setMethodS3("getForeground", "QuantArrayData", function(this, which=c("mean")) {
  which <- match.arg(which);
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(channels))
    channels <- 1:2;

  if (which == "mean") {
    fields <- c("ch1 Intensity", "ch2 Intensity")
  }

  # Assert that the fields do really exist.
  if (!all(fields %in% getFields(this))) {
    throw("The background estimates ", paste(fields, collapse=" and "), 
          " is not part of this ", data.class(this), " object.");
  }

  RGData(R=this[[fields[1]]], G=this[[fields[2]]], layout=getLayout(this));
})



setMethodS3("getForegroundSD", "QuantArrayData", function(this, slides=NULL, channels=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(channels))
    channels <- 1:2;

  res <- list();
  for (ch in channels) {
    field <- sprintf("ch%d Intensity Std Dev", as.integer(ch));
    res[[ch]] <- as.matrix(this[[field]][,slides]);
  }
  names(res) <- paste("fg std. dev. ch ", channels, sep="");
  
  res;
})


setMethodS3("getForegroundSNR", "QuantArrayData", function(this, slides=NULL, channels=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(channels))
    channels <- 1:2;
  
  fg <- list();
  for (ch in channels) {
    field <- sprintf("ch%d Intensity", as.integer(ch));
    fg[[ch]] <- as.matrix(this[[field]][,slides]);
  }
  
  fgSd <- getForegroundSD(this, slides=slides, channels=channels);
  
  res <- list();
  for (kk in seq(fg)) {
    res[[kk]] <- fg[[kk]] / fgSd[[kk]];
  }
  names(res) <- paste("fg SNR ch ", channels, sep="");
  res;
})

setMethodS3("getBackground", "QuantArrayData", function(this, which=c("mean")) {
  which <- match.arg(which);

  if (which == "mean") {
    fields <- c("ch1 Background", "ch2 Background")
  }

  # Assert that the fields do really exist.
  if (!all(fields %in% getFields(this))) {
    throw("The background estimates ", paste(fields, collapse=" and "), 
          " is not part of this ", data.class(this), " object.");
  }

  RGData(R=this[[fields[1]]], G=this[[fields[2]]], layout=getLayout(this));
})



setMethodS3("getBackgroundSD", "QuantArrayData", function(this, slides=NULL, channels=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(channels))
    channels <- 1:2;

  res <- list();
  for (ch in channels) {
    field <- sprintf("ch%d Background Std Dev", as.integer(ch));
    res[[ch]] <- as.matrix(this[[field]][,slides]);
  }
  names(res) <- paste("bg std. dev. ch ", channels, sep="");
  
  res;
})




setMethodS3("setLayout", "QuantArrayData", function(this, layout) {
  warning("For a QuantArrayData object it should not be necessary to set the layout explicitly. The layout is automatically calculated when the data is loaded from file.");
  NextMethod("setLayout", this);
})


# # setMethodS3("getArea", "QuantArrayData", function(this, slides=NULL, include=NULL, ...) {
# #   if (this$version <= 2)
# #     throw("getArea() is only supported for ", data.class(this), " version 3 and higher.");
# #   
# #   if (is.null(slides)) slides <- seq(nbrOfSlides(this));
# #   if (is.null(include)) include <- seq(nbrOfSpots(this));
# #   (this[["ch1 Area"]][include,slides] +
# #    this[["ch2 Area"]][include,slides]) / 2;
# # })


setMethodS3("anonymize", "QuantArrayData", function(this, ...) {
  # Anonymize a copy of the Layout object; others might use this one.
  layout <- clone(getLayout(this));
  anonymize(layout, ...);
  # Assume the same layout on all slides!
  slides <- 1:nbrOfSlides(this);
  this[["Name"]][,slides] <- getID(layout);
  setLayout(this, layout);
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
setMethodS3("getSpotPosition", "QuantArrayData", function(this, slides=NULL, index=NULL) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (any(index < 1) || any(index > nbrOfSpots(this)))
    throw("Argument 'index' is out of range.");
  if (is.null(index))
    index <- 1:nbrOfSpots(this);

  xField <- "X Location";
  yField <- "Y Location";
  if (!hasField(this, xField) || !hasField(this, yField)) {
    throw("This ", data.class(this), " object is missing the fields '", xField, "' and/or '", yField, "', which specifies the physical position of the spots.");
  }

  x <- this[[xField]][index,slides];
  y <- this[[yField]][index,slides];
  SpotPosition(x=x, y=y);
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
#   \item{yaxt, ...}{Parameters as accepted by \code{plot}.}
# }
#
# \value{Returns nothing.}
#
# @author
#
# \examples{
#  \dontrun{
#   # At the moment there is no QuantArray sample files in the package...
#   # qa <- QuantArrayData$read("quantarray123.txt", path=system.file("data-ex", package="aroma"))
#   # ...will use a GenePix sample file instead.
#   qa <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#
#   subplots(2)
#
#   opar <- par(bg="black")
#   plotSpatial(qa)
#   par(opar)
#
#   opar <- par(bg="black")
#   plotSpatial(qa, palette="blueyellow")
#   par(opar)
#  }
# }
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("plotSpatial", "QuantArrayData", function(this, what=NULL, slide=1, pch=20, yaxt=NULL, col="auto", palette="redgreen", A.range=c(0,16), M.range=c(-1,1), style=c("real", "classic"), ...) {
  style = match.arg(style);
  if (style == "classic")
    return(plotSpatial.MicroarrayData(this, what=what, slide=slide, col=col, ...));
  
  if (length(slide) != 1 || !is.numeric(slide))
    throw("Argument 'slide' must be an integer.");
  
  if (!is.null(what) && !is.element(what, getFieldNames(this)))
    throw("Argument 'what' is refering to an unknown field: ", what);
  
  if (is.null(col) || col == "auto") {
    if (is.null(what)) {
      raw <- getRawData(this, slides=slide);
      col <- MicroarrayData$createColors(raw$R,raw$G, type="RG", palette=palette, A.range=A.range, M.range=M.range);
    } else {
      col <- getColors(this, what, slide=slide);
    }
  }
  
  xy <- getSpotPosition(this, slides=slide);
  # Upper *left* corner is (0,0)
  x <- xy$x; y <- -xy$y; # Will also be used as xlab and ylab.
  
  if (missing(yaxt))
    yaxt0 <- "n";

  plot(x,y, pch=pch, yaxt=yaxt0, col=col, ...);
  
  if (missing(yaxt)) {
    # Add the y-axis with the correct ticks
    yaxp <- par("yaxp");
    yaxis <- seq(yaxp[1],yaxp[2], length=yaxp[3]+1);
    axis(2, at=yaxis, labels=-yaxis);
  }
})


setMethodS3("plotSpatial3d", "QuantArrayData", function(this, field=NULL, ...) {
  NextMethod("plotSpatial3d", this, field=field, ...);
})


setMethodS3("getArea", "QuantArrayData", function(this, slides=NULL, include=NULL, method=c("area", "diameter"), ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  fields <- getFieldNames(qa);

  method <- match.arg(method);
  found <- FALSE;
  if ((!found || method == "area") && is.element("ch1 Area", fields)) {
    a1 <- this[["ch1 Area"]][include,slides];
    a2 <- this[["ch2 Area"]][include,slides];
    found <- TRUE;
  }

  if ((!found || method == "diameter") && is.element("ch1 Diameter", fields)) {
    a1 <- this[["ch1 Diameter"]][include,slides];
    a2 <- this[["ch2 Diameter"]][include,slides];
    a1 <- pi*(a1/2)^2;
    a2 <- pi*(a2/2)^2;
    found <- TRUE;
  }

  if (!found) {
    a1 <- matrix(NA, nrow=length(include), ncol=length(slides));
    a2 <- NA;
  }
  
  as.matrix((a1 + a2) / 2);
})


############################################################################
# HISTORY:
# 2006-02-11
# o BUG FIX: Since (at least) R v2.1.1, the R internal readTableHead() takes
#   six arguments and not five.  This gave "Error in method(static, ...) : 
#   5 arguments passed to 'readTableHead' which requires 6".
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
# 2004-05-02
# o Added getForegroundSD(), getForegroundSNR(), getBackground(), and
#   getBackgroundSD().
# o Added Rdoc comments for the QA fields.
# o BUG FIX: Private read.table() generated an error in R v1.9.0.
# 2004-03-17
# o BUG FIX: The check for if *any* of the read unicoded characters
#   contained a hi-byte != 0 was only done on the first byte and not on
#   all(), which resulted in a warning "the condition has length > 1 and
#   only the first element will be used...".
# 2004-03-15
# o BUG FIX: Corrected a bug in as.RawData() complaining about 'this', which
#   was introduced in previous version.
# 2004-03-09
# o BUG FIX: write() was mistakenly accepting more than once slide at the
#   time. It was a typo and now an error is thrown is more slides are given.
# 2004-02-25
# o BUG FIX: Replaced broken isUnicoded@endian with attr(isUnicoded,
#   "endian"). Apparently, there's been a change in R making the former
#   invalid.
# 2003-10-26
# o Added support for *.gz files.
# 2003-09-10
# o Made the help text about read() that talks about Unicode support and how
#   it is done available again. It for mistakenly documented in a private
#   method.
# 2003-08-29
# o Added match.arg() where applicable.
# 2003-04-12
# o Added getBackground() and getForeground().
# 2002-12-21
# o Added plotSpatial3d().
# o When reading slides from files each slide is now named as the filename.
# 2002-12-05
# o BUG FIX: The get*() methods did not return a matrix if there was only
#   one slide.
# 2002-09-24
# o Changed the attribute 'path="."' to 'path=""' in read().
# o Update the Rdoc's for as.RawData().
# 2002-09-21
# o read() can now read one or several files specified by names or by a
#   pattern. This is identical to readAll(), which is now just calling
#   read() for backward compatibilities.
# 2002-08-22
# o Added support to read unicoded QuantArray files. This is done, by first
#   check if the file is unicoded and if it is the it checks whether it is
#   in Bigendian or Lowendian. Then the file is translated into a temporary
#   ASCII file by just excluding the highbyte and keeping the lowbytes.
#   This ASCII file is then read.
# o Added plotSpatial().
# o Added getSpotPosition().
# o Added getArea().
# 2002-08-21
# o Update readData() to not read integer's at all, but only assume
#   double's. Added a few comments about the Unicode format, which is
#   currently not supported by [R] nor com.braju.sma.
# 2002-08-18
# o Implemented readAll().
# o Since the 'Measurements' section in QuantArray files seems to contain
#   rows with tailing TAB's (that just should be ignored) read.table() fails
#   to read them. read.table() is making use of scan() and scan() has the
#   argument 'flush' which flushes such trailing cells, but it is not used
#   by read.table(). For this reason I created the static read.table()
#   method of QuantArrrayData which has the 'flush' argument.
# o Can read QuantArray v2 & v3 files.
# 2002-08-13
# o Making use of get AbsolutePath() in MicroarrayData.
# o Experience weird problems with the R CMD check. Maybe it is a perl bug?!
#   When removing the file it works find, but when including it *another*
#   file is reported to have syntax errors.
# 2002-08-12
# o Created!
############################################################################
