#########################################################################/**
# @RdocClass GenePixData
#
# @title "The GenePixData class"
#
# @synopsis
#
# \description{
#  @classhierarchy
# }
#
# \arguments{
#   \item{layout}{A @see "Layout" object specifying the spot layout of the
#    slides in this data set.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   For information about the fields in this class see Axon's
#   specifications of the GenePix File Formats at
#   \url{http://www.axon.com/GN_GenePix_File_Formats.html}.\cr
# }
#
# \section{Flagging feature-indicators}{
#   GenePix Pro can \emph{flags} individual spots according to
#   different criterias. The different flag statuses are
#   \emph{Good} (has value 100 in the GPR structure),
#   \emph{Bad} (-100), \emph{absent} (-75), and
#   \emph{Not Found} (-50). Unflagged spots have value 0.
#   A spot can only have one of these flags set at each time.
#
#   The \emph{Absent} flag is set \emph{automatically} if
#   there is a blank/empty entry in the GAL file (commonly
#   reflected in the @see "Layout" object associated with this
#   @see "GenePixData" object).
#
#   The \emph{Not Found} flag is set \emph{automatically} if
#   the spot/block alignment failed to find the spot. The 
#   position and diameter reported for such a spot are taken
#   from the initial values obtained from the rough gridding.
#   Note that the spot signals are still estimated.
#   It may happen that the spot is dislocated, but that the
#   user can visually pick it up and recenter the spot position
#   and ask GenePix Pro to resegment the spot. The result
#   may then be that GenePix flags the spot as \emph{Good}
#   or \emph{Not Found} if it still finds that the spot is too
#   weak for segmentation.
#
#   The \emph{Bad} and \emph{Good} flags are set
#   \emph{manually} by the user.
#   A typical scenario is that GenePix identified the spots
#   to have strong signals and segmented them well, but when
#   the user looks at the image she or he sees that two spots
#   are overlapping and have probably mixed their contents.
#   Then the user flags that spot as \emph{Bad}.
#   It can also be that there is a spot of interest, say a
#   negative control, and GenePix flags it as \emph{Not Found},
#   but the users looks at it and concludes that it is indeed
#   a high quality negative control spot and therefore flags
#   it as \emph{Good}.
#
#   Spots that have been flagged \emph{Absent}, \emph{Bad},
#   and \emph{Good} (not sure about this last one), will
#   \emph{not} be relabel by GenePix. Thus, GenePix will
#   only flag or unflag spots that are \emph{Not Found} or
#   unflagged.
#  
#   Note that any flag may be set or unset manually after the
#   automatic flagging has been done.
# }
#
# @author
#
# \references{
#   GenePix File Formats, 
#   \url{http://www.axon.com/GN_GenePix_File_Formats.html} and
#   \url{http://meetings.cshl.org/tgac/tgac/microarray.html}
# }
#
# \examples{
#   gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#
#   # Get the foreground and the background (and the layout)
#   raw <- getRawData(gpr)
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
#   The laser beam with wavelength 635nm is the red laser and excites the
#   Cy5 dye, and the one with wavelength 532nm is the green laser which
#   excites the Cy3 dye.
# }
#*/#########################################################################
setConstructorS3("GenePixData", function(layout=NULL, ...) {
  extend(MicroarrayData(layout=layout), "GenePixData",
    version=matrix()
  )
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#    Header fields
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getHeaderFields", "GenePixData", function(this, ...) {
  fields <- c();
  for (name in names(this)) {
    field <- this[[name]];
    if (!is.null(attr(field, "GPRHeaderField")) 
         && attr(field, "GPRHeaderField") == TRUE) {
      fields <- c(fields, name);
    }
  }
  fields;
}, protected=TRUE)


setMethodS3("getHeaderField", "GenePixData", function(this, pattern, ...) {
  headerFields <- getHeaderFields(this);
  if (is.null(headerFields))
    return(NULL);
  match <- regexpr(pattern, headerFields);
  idx <- which(match != -1);
  if (length(idx) == 0)
    return(NULL);
  if (length(idx) > 1)
    throw("Not a unique GPR header field name: ", pattern);
  name <- headerFields[idx];
  this[[name]];
}, protected=TRUE)


setMethodS3("setHeaderField", "GenePixData", function(this, name, value, ...) {
  tmp <- as.matrix(value);
  attr(tmp, "GPRHeaderField") <- TRUE;
  this[[name]] <- tmp;
}, protected=TRUE)





setMethodS3("append", "GenePixData", function(this, other, ...) {
  # Turn off virtual fields so getCircularity() etc is not called.
  attr(this, "disableGetMethods") <- TRUE;
  attr(other, "disableGetMethods") <- TRUE;

  thisFields <- getFields(this);
  otherFields <- getFields(other);

  # i) Add fields to other that exists in this but not other...
  missingInOther <- setdiff(thisFields, otherFields);
  for (field in missingInOther) {
    value <- this[[field]];
    value[] <- NA;
    other[[field]] <- value;
  }
  # ii) ...and vice versa
  missingInThis <- setdiff(otherFields, thisFields);
  for (field in missingInThis) {
    value <- other[[field]];
    value[] <- NA;
    this[[field]] <- value;
  }

  this$.fieldNames <- unique(c(getFieldNames(this), getFieldNames(other)));

  # Append all fields names
  NextMethod("append");

  # Append all GPRHeaderField:s by assuming they are all
  # one dimensional, i.e. column vectors.
  for (name in getHeaderFields(this)) {
    a <- as.matrix(this[[name]]);
    b <- as.matrix(other[[name]]);
    if (!identical(nrow(a), nrow(b))) {
      t <- rep(NA, length=nrow(a));
      t[1:nrow(b)] <- b;
      b <- t;
    }
    setHeaderField(this, name, cbind(a,b));
  }

  # Append the version fields too.
  if (inherits(other, "GenePixData"))
    this$version <- cbind(this$version, other$version)
  else 
    this$version <- cbind(this$version, NA);

  invisible(this);
})


setMethodS3("getKnownHeadersWithType", "GenePixData", function(static, settings=NULL, verbose=FALSE, ...) {
  # This methods returns a list of all known headers and their known
  # R datatype to be used when read by read.table() etc. 

  # NOTE: Some datatypes are currently unknown and are therefore set to
  # the read generic type NA, which is always understood. Hence, for
  # writing some of these types has to be updated before writing.

  # Regular expressions of headers. Watch out for "."'s and "+"'s!
  # Note that the "Log Ratio" field sometimes contains "Error", i.e. we
  # can't use "double".

  channelNames <- getChannelNames(static, settings=settings);
  incl <- !is.na(channelNames);
  channels <- rep(NA, length=4);
  channels[incl] <- channelNames[incl];

  if (identical(verbose, TRUE)) {
    cat("Identified channel names: ", paste(channels, collapse=", "), "\n", sep="");
  }

  # GenePix Pro x.x.x.x (GenePix Results Format Version 1.0?)
  # Source: http://www.axon.com/gn_GPR_Format_History.html
  knownHeaders <- c(
    "Block"="integer", 
    "Column"="integer",
    "Row"="integer",
    "Name"="character",
    "ID"="character",
    "X"="integer",
    "Y"="integer",
    "Dia\\."="integer",
    "F<channel[1]> Median"="integer",
    "F<channel[1]> Mean"="integer",
    "F<channel[1]> SD"="integer",
    "B<channel[1]> Median"="integer",
    "B<channel[1]> Mean"="integer",
    "B<channel[1]> SD"="integer",
    "% > B<channel[1]>\\+1SD"="integer",
    "% > B<channel[1]>\\+2SD"="integer",
    "F<channel[1]> % Sat\\."="integer",
    "F<channel[2]> Median"="integer",
    "F<channel[2]> Mean"="integer",
    "F<channel[2]> SD"="integer",
    "B<channel[2]> Median"="integer",
    "B<channel[2]> Mean"="integer",
    "B<channel[2]> SD"="integer",
    "% > B<channel[2]>\\+1SD"="integer",
    "% > B<channel[2]>\\+2SD"="integer",
    "F<channel[2]> % Sat\\."="integer",
    "Ratio of Medians"="double",
    "Ratio of Means"="double",
    "Median of Ratios"="double",
    "Mean of Ratios"="double",
    "Ratios SD"="double",
    "Rgn Ratio"="double",
    "Rgn R"="double",
    "F Pixels"="integer",
    "B Pixels"="integer",
    "Sum of Medians"="integer",
    "Sum of Means"="integer",
    "Log Ratio"=NA,
    "F<channel[1]> Median - B<channel[1]>"="integer",
    "F<channel[2]> Median - B<channel[2]>"="integer",
    "F<channel[1]> Mean - B<channel[1]>"="integer",
    "F<channel[2]> Mean - B<channel[2]>"="integer",
    "Flags"="integer",
    "Normalize"="integer",
    "F<channel[3]> Median - B<channel[3]>"="integer",
    "F<channel[4]> Median - B<channel[4]>"="integer",
    "F<channel[3]> Mean - B<channel[3]>"="integer",
    "F<channel[4]> Mean - B<channel[4]>"="integer",
    "SNR <channel[3]>"="double",
    "F<channel[3]> Total Intensity"="integer",
    "Index"="integer",
    "User Defined"=NA
  );

  # GenePix Pro 3.0.6.x (GenePix Results Format Version 1.4)
  # Source: http://www.axon.com/gn_GPR_Format_History.html
  knownHeaders <- c(knownHeaders, 
    # I am not sure about when four-colour fields was first introduced
    # but I add it here. /HB 2003-08-29
    "F<channel[3]> Median"="integer", 
    "F<channel[3]> Mean"="integer", 
    "F<channel[3]> SD"="integer", 
    "B<channel[3]> Median"="integer",
    "B<channel[3]> Mean"="integer",
    "B<channel[3]> SD"="integer",
    "% > B<channel[3]>\\+1SD"="integer",
    "% > B<channel[3]>\\+2SD"="integer",
    "F<channel[3]> % Sat\\."="integer",
    "F<channel[4]> Median"="integer", 
    "F<channel[4]> Mean"="integer",
    "F<channel[4]> SD"="integer",
    "B<channel[4]> Median"="integer",
    "B<channel[4]> Mean"="integer",
    "B<channel[4]> SD"="integer",
    "% > B<channel[4]>\\+1SD"="integer",
    "% > B<channel[4]>\\+2SD"="integer",
    "F<channel[4]> % Sat\\."="integer",
    "F<channel[3]> Median - B<channel[3]>"="integer",
    "F<channel[4]> Median - B<channel[4]>"="integer",
    "F<channel[3]> Mean - B<channel[3]>"="integer",
    "F<channel[4]> Mean - B<channel[4]>"="integer"
  );

  # GenePix Pro 4.0.1.x (GenePix Results Format Version 2.0)
  # Source: http://www.axon.com/gn_GPR_Format_History.html
  knownHeaders <- c(knownHeaders, 
    "Ratio of Medians \\(<channel[1]>/<channel[2]>\\)"="double", # All ratio columns were re- 
    "Ratio of Medians \\(Ratio/2\\)"="double", # named on a per ratio basis,
    "Ratio of Medians \\(Ratio/3\\)"="double", # e.g., "Ratio of Medians" in
    "Ratio of Means \\(<channel[1]>/<channel[2]>\\)"="double",   # GenePixPro 3.0.6 is changed
    "Ratio of Means \\(Ratio/2\\)"="double",   # to "Ratio of Medians       
    "Ratio of Means \\(Ratio/3\\)"="double",   # (<channel[1]>/<channel[2]>)". The calculation
    "Median of Ratios \\(<channel[1]>/<channel[2]>\\)"="double", # did not change.            
    "Median of Ratios \\(Ratio/2\\)"="double",
    "Median of Ratios \\(Ratio/3\\)"="double",
    "Mean of Ratios \\(<channel[1]>/<channel[2]>\\)"="double",  
    "Mean of Ratios \\(Ratio/2\\)"="double",
    "Mean of Ratios \\(Ratio/3\\)"="double",
    "Ratios SD \\(<channel[1]>/<channel[2]>\\)"="double",       
    "Ratios SD \\(Ratio/2\\)"="double",
    "Ratios SD \\(Ratio/3\\)"="double",
    "Rgn Ratio \\(<channel[1]>/<channel[2]>\\)"="double",       
    "Rgn Ratio \\(Ratio/2\\)"="double",
    "Rgn Ratio \\(Ratio/3\\)"="double",
    "Rgn R. \\(<channel[1]>/<channel[2]>\\)"="double", 	      
    "Rgn R. \\(Ratio/2\\)"="double",
    "Rgn R. \\(Ratio/3\\)"="double",
    "Log Ratio \\(<channel[1]>/<channel[2]>\\)"=NA,
    "Log Ratio \\(Ratio/2\\)"=NA,
    "Log Ratio \\(Ratio/3\\)"=NA,
    "Normalize"="integer"
  );

  # GenePix Pro 4.1.1.x (GenePix Results Format Version 3.0)
  # Source: http://www.axon.com/gn_GPR_Format_History.html
  knownHeaders <- c(knownHeaders, 
    "F<channel[1]> Total Intensity"="integer", 
    "F<channel[2]> Total Intensity"="integer", 
    "SNR <channel[1]>"="double",
    "SNR <channel[2]>"="double"
  );

  # GenePix Pro 5.0.0.x (GenePix Results Format Version 3.0)
  # Source: http://www.axon.com/gn_GPR_Format_History.html
  knownHeaders <- c(knownHeaders, 
    "Negative Control"=NA,
    "B<channel[1]>"=NA,             
    "B<channel[2]>"=NA,            
    "F<channel[2]> CV"=NA,
    "F<channel[1]> CV"=NA,
    "B<channel[2]> CV"=NA,
    "B<channel[1]> CV"=NA,
    "Circularity"="integer",
    "Autoflag"=NA
  );

  # Give all "channel" fields their correct names
  names(knownHeaders) <- convertFieldNames(static, names(knownHeaders), channelNames=channels);

  if (identical(verbose, TRUE)) {
    cat("Header considered to be known: ", paste(names(knownHeaders), collapse=", "), "\n", sep="");
  }

  knownHeaders;
}, private=TRUE, static=TRUE);



#########################################################################/**
# @RdocMethod readHeader
#
# @title "Reads a GenePix Results (GPR) file header"
#
# \description{
#   @get "title". This method will parse the header of a GPR file and
#   extract all header information including the column names of the
#   preceeding data section.
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
# \value{Returns a @list of the header fields, but also the names and
#   the "guessed" classes of the columns to be read.}
#
# @author
#
# \details{
#   This method is used internally by the @seemethod "read" method.
#
#   There are many difference version of the GenePix image analysis 
#   software and almost equally number of different versions of the
#   GenePix Results file format. All example data available at
#   \url{http://www.axon.com/gn_GPR_Format_History.html} have been
#   tested and can successfully be read using this method.
# }
#
# \references{
#   GenePix File Formats:
#   \url{http://www.axon.com/GN_GenePix_File_Formats.html},
#   \url{http://www.axon.com/gn_GPR_Format_History.html}
# }
#
# \examples{
#   header <- GenePixData$readHeader("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#   str(header)
# }
#
# \seealso{
#   To read a GenePix Results file including the data section 
#   see @seemethod "read".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("readHeader", "GenePixData", function(static, filename, path=NULL, verbose=FALSE, ...) {
  filename <- Arguments$getReadablePathname(filename, path);  

  if (verbose) cat("Reading the header of file ", filename, "...", sep="");

  # Support gzip'ed files too.
  if (regexpr("[.]gz$", filename) != -1) {
    tmpname <- tempfile();
    n <- gunzip(filename, tmpname);
    filename <- tmpname;
    on.exit({
      file.remove(tmpname);
    });
  }
  
  fh <- file(filename, "r"); 
  on.exit({
    if (inherits(fh, "connection") && isOpen(fh)) {
      close(fh);
    }
  }, add=TRUE);


  #=========================================================================
  # GPR Header
  #=========================================================================

  #-------------------------------------------------------------------------
  # Line 1: "ATF   1.0"
  # File type and version number.
  #-------------------------------------------------------------------------
  lines <- scan(fh, what=list("",""), sep="\t", quote="", nlines=1, strip.white=TRUE, flush=TRUE, quiet=TRUE);
  fileFmt <- lines[[1]];
  if (fileFmt != "ATF") {
    warning(paste("Unknown file format: ", fileFmt, "\n", sep=""));
  }

  # Currently the ATF version number can be "1.0" or "1".
  version <- lines[[2]];
  if (version != "1.0" && version != "1") {
    warning(paste("Unknown ATF file format version: ", fileFmt, " ", version, "\n", sep=""));
  }

  if (verbose)
    cat(fileFmt, version, "\n");

  #-------------------------------------------------------------------------
  # Line 2: "24   43"
  # Number of optional header records and number of data fields (columns).
  #-------------------------------------------------------------------------
  lines <- scan(fh, what=list(0,0), sep="\t", quote="", nlines=1, strip.white=TRUE, flush=TRUE, quiet=TRUE);
  nbrOfOptionalHeaders <- lines[[1]];
  nbrOfColumns <- lines[[2]];

  if (verbose)
    cat(nbrOfOptionalHeaders, nbrOfColumns, "\n");

  if (verbose)
    cat("Reading", nbrOfOptionalHeaders, "optional header records:\n");

  settings <- c();
  for (k in seq(nbrOfOptionalHeaders)) {
    #-----------------------------------------------------------------------
    # Line 3:2+nbrOfOptionalHeaders: 
    # "FieldName=Field Value"
    #-----------------------------------------------------------------------
    lines <- scan(fh, what=list(""), sep="\t", quote="\"", nlines=1, strip.white=TRUE, flush=TRUE, quiet=TRUE);
    record <- strsplit(lines[[1]], "=");
    field.name <- record[[1]][1];
    field.value <- record[[1]][2];
    settings[field.name] <- field.value;
    if (verbose)
      cat(field.name, "=", field.value, "\n", sep="");
  }

  #-------------------------------------------------------------------------
  # Line 2+nbrOfOptionalHeaders+1:end
  #-------------------------------------------------------------------------
  # Regular expressions of headers. Watch out for "."'s and "+"'s!
  knownHeaders <- getKnownHeadersWithType(static, settings=settings, verbose=verbose);
  # Prescan the headers to identify version
  if (verbose)
    cat("Reading field names...");
  # scan() can't read column names that contains non-ASCII 0-127 characters!
  header <- readLines(con=fh, n=1);
  close(fh); fh <- NULL;
  
  # Remove leading blanks
#  header <- gsub("^[ ]*", "", header);  # Does not work. Why? /HB 020807
  # Remove trailing blanks
  header <- gsub("[ ]*$", "", header);
  # Remove quotation marks
  header <- gsub("\"", "", header);

  # Extract the column names
  header <- unlist(strsplit(header, "\t"));
  if (verbose)
    cat("ok\n");

  unknownHeaders <- c();
  colClasses <- c();
  for (k in seq(along=header)) {
    colNames <- names(knownHeaders);
    typeIdx <- which(sapply(seq(colNames), FUN=function(i) regexpr(colNames[i], header[k]) == 1));
    if (length(typeIdx) == 0) {
      unknownHeaders <- c(unknownHeaders, k);
      colClasses <- c(colClasses, NA);
    } else {
      typeIdx <- typeIdx[1];  # Bullet proof!
      colClasses <- c(colClasses, knownHeaders[typeIdx]);
    }
  }

  if (length(unknownHeaders) > 0) {
    unknown <- paste(header[unknownHeaders], collapse=", ");
    warning(paste(sep="", "GenePix file format warning: Some of the field names are not recognized. It might be because it is a new/old version. The unknown fields are: ", unknown, "."));
  }

  gprHeader <- list(
    version              = as.matrix(version),
    colNames             = header,
    colClasses           = colClasses,
    nbrOfOptionalHeaders = nbrOfOptionalHeaders
  );

  # Make all "settings" into GPR fields like the other data fields.
  for (name in names(settings)) {
    tmp <- settings[[name]];
    tmp <- strsplit(tmp, split="\t")[[1]];
    tmp <- as.matrix(tmp);
    attr(tmp, "GPRHeaderField") <- TRUE;
    gprHeader[[name]] <- tmp;
  }

  if (verbose) cat("ok\n", sep="");
  
  gc(); # To minimize memory usage!

  gprHeader;
}, static=TRUE);  # readHeader()



#########################################################################/**
# @RdocMethod read
#
# @title "Reads a GenePix Results Data file"
#
# \description{
#   Reads the GenePix file format GPR (GenePix Results Format).
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
# \value{Returns a @see "GenePixData" object.}
#
# @author
#
# \details{
#   There are many difference version of the GenePix image analysis 
#   software and almost equally number of different versions of the
#   GenePix Results file format. All example data available at
#   \url{http://www.axon.com/gn_GPR_Format_History.html} have been
#   tested and can successfully be read using this method.
# }
#
# \references{
#   GenePix File Formats:
#   \url{http://www.axon.com/GN_GenePix_File_Formats.html},
#   \url{http://www.axon.com/gn_GPR_Format_History.html}
# }
#
# \examples{
#   gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
# }
#
# \seealso{
#   To write a slide to a GenePix Results file see @seemethod "write".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("readOneFile", "GenePixData", function(this, filename, path=NULL, solve=FALSE, verbose=FALSE, ...) {
  filename <- Arguments$getReadablePathname(filename, path);  

  if (verbose) cat("Reading file ", filename, "...", sep="");

  # Support gzip'ed files too.
  if (regexpr("[.]gz$", filename) != -1) {
    tmpname <- tempfile();
    n <- gunzip(filename, tmpname);
    filename <- tmpname;
    on.exit(file.remove(tmpname));
  }
  
  #=========================================================================
  # GPR Header
  #=========================================================================
  gprHeader <- readHeader(this, filename=filename, path=NULL, verbose=verbose);
  nbrOfOptionalHeaders <- gprHeader$nbrOfOptionalHeaders;
  colNames   <- gprHeader$colNames;
  colClasses <- gprHeader$colClasses;


  #=========================================================================
  # GPR Data
  #=========================================================================
  # Note that the "Log Ratio" field sometimes contains "Error", i.e. we
  # can't use "double".
  if (verbose)
    cat("Reading result table...");
  # NOTE: read.table(fh,...) does NOT work! Probably a bug in the skip part!
  # /Henrik 2001-03-07

  # Try different types of quotings
  df <- NULL;
  quotes <- c("", "\"", "'\"", "'");
  for (quote in quotes) {
    tryCatch({
      if (solve == TRUE) {
#        if (verbose)
#          cat("Number of columns expected: ", length(colClasses), ".\n", sep="");
        # It is not possible to set colClasses <- NULL, which gives an error.
        df <- read.table(filename, sep="\t", quote=quote, comment.char="", 
           skip=nbrOfOptionalHeaders+3, header=FALSE, na=c("NA", "NaN", "Error"));
      } else {
        df <- read.table(filename, sep="\t", quote=quote, comment.char="", 
           skip=nbrOfOptionalHeaders+3, colClasses=colClasses, header=FALSE, 
                                                        na=c("NA", "NaN", "Error"));
      }
    }, error = function(ex) {
      if (quote == quotes[length(quotes)]) {
        if (verbose)
          cat("failed!\n");
        stop("Failed to read GPR result table: ", ex$message);
      }
    })
    if (!is.null(df))
      break;
  }
  if (verbose)
    cat(nrow(df), "lines. ok!\n");
  colnames(df) <- colNames;

  # "Name" and "ID" becomes factors, which are "harder" to work with.
  # By changing them to vectors the size increases (just) a little bit.
  # Change all factor levels that should be character's.
  for (k in which(colClasses == "character"))
    df[[k]] <- as.character(df[[k]]);

  # Change all factor levels that should be double's.
  for (k in which(colClasses == "double"))
    df[[k]] <- as.numeric(as.character(df[[k]]));
  
  if (verbose)
    cat("Successfully read file.\n");

  #-------------------------------------------------------------------------
  # Extracting additional field to make life easier...
  #-------------------------------------------------------------------------
  nbrOfGrids <- max(df$Block);
  nspot.r    <- max(df$Row);
  nspot.c    <- max(df$Column);

  #-------------------------------------------------------------------------
  # Figuring out the layout structure of the microarray
  #-------------------------------------------------------------------------
  # Tries to figure out the number of rows of grids from the x coordinates.

  # Number of spots in each grid
  gridSize <- nspot.r * nspot.c;

  # Number of genes
  nbrOfSpots <- nbrOfGrids * gridSize;
  # The previous line is a little bit safer than the following line
  # if some rows are missing:
  # nbrOfSpots <- length(df[,1]);

  #-------------------------------------------------------------------------
  # Make sure that the spots are ordered block -> spot row -> spot column
  #-------------------------------------------------------------------------
  # Calculate a hashcode based on (Block,Row,Column):
  observedSpots <- df[,c("Block", "Column", "Row")];
  # Range of integers: [-(2^31-1),+(2^31-1)] + {NA} > [-1000^3,+1000^3]
  ROW.BASE   <- as.integer(1000);
  BLOCK.BASE <- as.integer(1000*ROW.BASE);
  observedHash <- as.integer(BLOCK.BASE*observedSpots$Block + ROW.BASE*observedSpots$Row + observedSpots$Column);
  rm(observedSpots);
  order <- order(observedHash);
  if (!identical(1:length(order), order)) {
    # Not in order, reorder...
    msg <- sprintf("The spots are not listed in a ascending (Block,Row,Column) order. Will reorder them accordingly. Please be aware that spatial plots and spatial analysis that does not rely on the true physical positions, but assumed position according to GenePix's (Block,Row,Column) index may not be correct.\n");
    warning(msg);
    if (verbose) cat(msg, "\n");
    df <- df[order,];
    observedHash <- observedHash[order];
  }

  #-------------------------------------------------------------------------
  # If "missing" spots exist, insert NAs for those...
  #-------------------------------------------------------------------------
  if (nbrOfSpots != nrow(df)) {
    nbrOfMissingSpots <- nbrOfSpots - nrow(df);

    msg <- sprintf("Printtip groups (aka blocks & grids) were not complete. Will reconstruct them by inserting in total %d \"missing\" spots.", as.integer(nbrOfMissingSpots));
    warning(msg);

    if (verbose) cat(msg, "\n");

    # Some rows in the data file is missing. Tries to reconstruct those.
    # 1. Identify missing rows by looking at "Block", "Column" and "Row".
    # a) we expect the following values
    expectedSpots <- data.frame(
      Block  = rep(1:nbrOfGrids, each=gridSize),
      Column = rep(1:nspot.c, times=nbrOfSpots %/% nspot.c),
      Row    = rep(rep(1:nspot.r, each=nspot.c), times=nbrOfGrids)
    );
    expectedHash <- BLOCK.BASE*expectedSpots$Block + ROW.BASE*expectedSpots$Row + expectedSpots$Column;

    missingRows <- which(!(expectedHash %in% observedHash));
    if (length(missingRows) != nbrOfMissingSpots)
      throw("Strange, the number of missing spots does not match the number of (Block,Row,Column) indices that was identified to be missing. Please report this error to the package author.");

    idx <- insert(1:nrow(df), missingRows);
    df <- df[idx,];
    df[missingRows, c("Block", "Column", "Row")] <- expectedSpots[missingRows, c("Block", "Column", "Row")];
    rm(expectedSpots, missingRows);

    if (verbose)
      cat("done\n");
  }

  #-------------------------------------------------------------------------
  # Figure out where the spots changes row
  #-------------------------------------------------------------------------
  # 1a. Indices of the *first* spot in each grid
  idx <- which(df$Row == 1 & df$Column == 1);

  for (offset in 0:(gridSize-1)) {
    # 1b. The x coordinate for all such spots
    x <- df$X[idx+offset];
  
    # 1c. The horizontal (x) distance to the next spot
    dx <- diff(c(x,0));
  
    # 1d. If the distance to the next spot is negative, we know that we
    #     have started over on a new row, i.e. changed row. The number of
    #     such changes (including the "fake" last row) is equal to number 
    #     of grid rows.
    ngrid.r <- sum(dx < 0);

    # 1e. If there were NAs (basically from inserted spots above), the sum
    #     is NA and we can not make use of this 'nspot.r' guess. Shift one
    #     spot within a grid and try again, otherwise finish.
    if (!is.na(ngrid.r))
      break;
  }
  rm(idx,x,dx,offset);

  # The number of grid columns.
  ngrid.c <- max(df$Block) / ngrid.r;

  # Assert that out calculation are correct.
  if (round(ngrid.c) != ngrid.c)
    warning("Number of grid rows seems to be wrong! Please report this error to henrikb@braju.com.\n");

  # Make sure the name's and the id's are not factors!
  name <- as.character(df[,"Name"]);
  id <- as.character(df[,"ID"]);

  # Trim names and ids (because names with trailing newlines have been seen)
  name <- trim(name);
  id <- trim(id);

  layout <- Layout(ngrid.r, ngrid.c, nspot.r, nspot.c, name=name, id=id)
  gpr <- GenePixData(layout=layout)
  gpr$version <- gprHeader$version;
  for (kk in seq(along=gprHeader)) {
    value <- gprHeader[[kk]];
    type <- attr(value, "GPRHeaderField");
    if (identical(type, TRUE)) {
      name <- names(gprHeader)[kk];
      gpr[[name]] <- value;
    }
  }
  
  gpr$.fieldNames <- names(df);
  for (field in names(df)) {
    gpr[[field]] <- as.matrix(df[[field]]);
    df[[field]] <- NULL; # Save memory
  }

  if (verbose) cat("ok\n", sep="");
  
  rm(df); gc(); # To minimize memory usage!

  gpr;
}, private=TRUE, static=TRUE);



#########################################################################/**
# @RdocMethod write
#
# @title "Write a GenePix Results Data file"
#
# \description{
#  Writes the GenePixData object to a file with the GenePix file format
#  GPR (GenePix Results Format).
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file to be written.}
#   \item{path}{The path to the file.}
#   \item{slide}{An integer specifying which slide to be written to file. 
#     Currently, only one slide at the time can be written.}
#   \item{overwrite}{If @TRUE, already an existing file will be overwritten. 
#       Otherwise an Exception will be thrown.}
#   \item{...}{Arguments passed to \code{write.table}.}
# }
#
# \value{Returns nothing.}
#
# @author
#
# \references{
#   GenePix File Formats, 
#   \url{http://www.axon.com/GN_GenePix_File_Formats.html}
# }
#
# \examples{
#   gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#
#   # Writes the GenePix Results Data to a file named "temp.gpr". Note
#   # that this file won't be exactly the same since a few lines,
#   # specifying for instance the creator of the file, will be added.
#   write(gpr, "temp.gpr", overwrite=TRUE)
#   file.show("temp.gpr")
# }
#
# \seealso{
#   To read one or more GenePix Results files
#   see @seemethod "read".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("write", "GenePixData", function(this, filename, path=NULL, slide=1, overwrite=FALSE, ..., verbose=FALSE) {
  filename <- Arguments$getWritablePathname(filename, path, mustNotExist=!overwrite);  

  slide <- validateArgumentSlide(this, slide=slide);

  if (!is.element("Creator", names(this))) {
    creator <- paste("com.braju.sma, http://www.braju.com/R/,", date());
    creator <- matrix(creator, nrow=1, ncol=nbrOfSlides(this));
    setHeaderField(this, "Creator", creator);
  }
  
  fh <- file(filename, "w");
  on.exit({
    # Will be called even on an error.
    if (inherits(fh, "connection") && isOpen(fh)) {
      close(fh);
    }
  }, add=TRUE);
  
  # Write ATF header
  cat("ATF\t", this$version[,slide], "\n", file=fh);
  cat(file=fh, length(getHeaderFields(this)), "\t", nbrOfFields(this), "\n", sep="");

  # Write GPR headers
  for (name in getHeaderFields(this)) {
    setting <- this[[name]][,slide];
    setting <- paste(setting, collapse="\t");
    setting <- paste("\"", name, "=", setting, "\"\n", sep="");
    cat(file=fh, setting, append=TRUE);
  }

  # Get known headers and their types
  settings <- c(Type=getHeaderField(this, "^Type.*"), Wavelengths=getHeaderField(this, "^Wavelengths.*"));
  knownHeaders <- getKnownHeadersWithType(this, settings=settings, verbose=verbose);
  
  # Get the column classes for each data field.
  colNames <- getFieldNames(this);
  knownNames <- names(knownHeaders);
  colClasses <- rep(NA, length(colNames));
  for (k in 1:length(colNames)) {
    found <- FALSE;
    for (l in seq(knownNames)) {
      found <- (regexpr(knownNames[l], colNames[k]) != -1);
      if (found) break;
    }
    if (found)
      colClasses[k] <- knownHeaders[l];
  }
  rm(knownNames, knownHeaders);
  
  stringCols <- (!is.na(colClasses) & colClasses == "character");

  # Write GPR data header
  colNames <- paste(paste("\"", colNames, "\"", sep=""), collapse="\t");
  cat(colNames, "\n", file=fh, append=TRUE, sep="");
  rm(colNames);
  
  # Create the separator strings
  sep <- "\t";
  eol <- "\n";
  sepStrings <- c(rep(sep, nbrOfFields(this)-1), eol);
  rm(sep, eol);

  # NOTE: To minimize memory usage, write chunks of 512 rows at the time
  chunkSize <- 512;
  nrows <- nbrOfSpots(this);
  nchunks <- nrows %/% chunkSize + 1;
  from <- 1;
  for (k in 1:nchunks) {
    gc();
    
    to <- min(from+chunkSize-1, nrows);
    if (verbose == TRUE)
      cat("Writing data rows ", from, ":", to, "\n", sep="");
#    cat(System$currentTimeMillis(),"\n")
    
    df <- extract(this, slide=slide, index=from:to);
  
    # Convert each column to its correct data type
    for (k in which(!is.na(colClasses[k]))) {
      if (colClasses[k] == "integer")
        df[,k] <- as.integer(df[,k])
      else if (colClasses[k] == "double")
        df[,k] <- as.double(df[,k])
      else if (colClasses[k] == "character")
        df[,k] <- as.character(df[,k])
    }
    
    # Remember which values are NA's
    nas <- is.na(df);
  
    df <- as.matrix(df);
    
    # Set the NA's to "Error"
    df[nas] <- "Error";
    
    # Strip blanks from non-strings; saves about 20% of the file size.
    df[,!stringCols] <- gsub(" ", "", df[,!stringCols]);

    df <- paste(c(t(df)), sepStrings, sep="", collapse="");
    writeLines(df, con=fh, sep="");
    
    from <- to + 1;
  } # for (k in 1:nchunks)
  
  # Calling write.table() with add a lot of unnecessary blanks to numeric columns.
  # Don't use it! /HB 020807
  #  write.table(df, file=fh, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", na="Error", append=TRUE, ...);
})


#########################################################################/**
# @RdocMethod read
#
# @title "Reads one or several GenePix files into one GenePixData object"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{A vector of filenames. Either \code{pattern} or \code{filename} must be specified.}
#   \item{path}{A string (or an optional vector of paths if \code{filename} is specified) to the files.}
#   \item{pattern}{A pattern string for matching filenames. Either \code{pattern} or \code{filename} must be specified.}
#   \item{verbose}{If @TRUE, information will printed out during
#                  the reading of the file.}
# }
#
# \value{Returns a @see "GenePixData" object.}
#
# @author
#
# \examples{
#   # Reads one GenePix data files
#   gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#
#   # Reads two GenePix data files
#   filenames <- c("gpr123.gpr", "gpr123.gpr");
#   gpr <- GenePixData$read(filenames, path=system.file("data-ex", package="aroma"))
#
#   # Reads all GenePix data files
#   pattern <- c(".*\\.gpr");
#   gpr <- GenePixData$read(pattern=pattern, path=system.file("data-ex", package="aroma"))
# }
#
# \seealso{
#   To write a slide to a GenePix Results file see @seemethod "write".
#   For pattern formats see @see "base::list.files".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("read", "GenePixData", function(this, filename=NULL, path=NULL, pattern=NULL, verbose=FALSE, ...) {
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
    tmp <- GenePixData$readOneFile(filename=filename[k], path=path[k], verbose=verbose, ...);
    slidename <- basename(filename[k]);
    setSlideNames(tmp, slidename);
    if (is.null(res))
      res <- tmp
    else {
      append(res, tmp); # Will assert same layout.
    }
  }

  gc();

  res;
}, static=TRUE);


# Kept for backward compatibility.
setMethodS3("readAll", "GenePixData", function(...) {
  GenePixData$read(...);
}, private=TRUE, static=TRUE, deprecated=TRUE)




setMethodS3("setLayout", "GenePixData", function(this, layout, ...) {
  warning("For a GenePixData object it should not be necessary to set the layout explicitly. The layout is automatically calculated when the data is loaded from file.");
  NextMethod("setLayout", this);
})



#########################################################################/**
# @RdocMethod getRawData
#
# @title "Gets the raw intensites from the GPR structure"
#
# \description{
#  Extracts the red and green spot intensitites (both signal and background)
#  from the GenePixData object and returns a @see "RawData" object.
#
#  \emph{Note: From com.braju.sma v0.59, this method returns the
#        *median* foreground and background estimates. Previous 
#        versions returned the *mean* of dito.}
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{Specifying which slides to be extracted. If @NULL, all slides are considered.}
#   \item{fg}{If \code{"mean"}, the mean foreground intensities are returned. If \code{"median"}, the median foreground intensities are returned.}
#   \item{bg}{If \code{"mean"}, the mean background intensities are returned. If \code{"median"}, the median background intensities are returned.}
#   \item{channels}{A @vector of length two specifying which two channels 
#     (wavelengths) to be extracted in case the data contains more than two
#      channels (colours).}
# }
#
# \value{Returns a @see "RawData" object containing the specified slides.}
#
# \details{
#   The R and Rb channels will come from the F635* and B635* fields, and
#   the G and Gb channels will come from the F532* and B532* fields.
#   To swap the channels just use dyeSwap().
#
#   From com.braju.sma v0.30 the \code{dye.swap} argument has been removed.
#   Any dye swapping should be done bye using the dyeSwap() method.
# }
#
# @author
#
# \references{
#   GenePix File Formats, 
#   \url{http://www.axon.com/GN_GenePix_File_Formats.html}
# }
#
# \examples{
#   gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#
#   # Gets the foreground and the background
#   raw <- getRawData(gpr)
# }
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getRawData", "GenePixData", function(this, slides=NULL, fg="median", bg="median", channels=1:2, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  channels <- validateArgumentChannels(this, channels=channels);

  #-------------------------------------------------------------------------
  # RGs as in the sma library.
  #-------------------------------------------------------------------------
  fg <- getForeground(this, which=fg, slides=slides, channels=channels);
  bg <- getBackground(this, which=bg, slides=slides, channels=channels);

  RawData(R=fg$R, G=fg$G, Rb=bg$R, Gb=bg$G, layout=getLayout(this), 
                                           extras=this$.extras)
});


setMethodS3("as.RawData", "GenePixData", function(this, ...) {
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
#   \item{yaxt, ...}{Parameters as accepted by \code{plot}.}
# }
#
# \value{Returns nothing.}
#
# @author
#
# \examples{
#   gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#
#   subplots(2)
#
#   opar <- par(bg="black")
#   plotSpatial(gpr)
#   par(opar)
#
#   opar <- par(bg="black")
#   plotSpatial(gpr, palette="blueyellow")
#   par(opar)
# }
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("plotSpatial", "GenePixData", function(this, what=NULL, slide=1, include=NULL, exclude=NULL, pch=20, yaxt=NULL, col="auto", palette="redgreen", A.range=c(0,16), M.range=c(-1,1), look=c("real", "classic"), cex=NULL, style=NULL, ...) {
  cmd <- NULL;
  if (!is.null(style) && is.element(style, c("highlight", "points", "text"))) {
    cmd <- style;
    lastPlot <- Device$getPlotParameters();
    what <- lastPlot$what;
    slide <- lastPlot$slide;
    look <- lastPlot$look;
  }

  look = match.arg(look);
  if (look == "classic")
    return(plotSpatial.MicroarrayData(this, what=what, slide=slide, col=col, ...));
  
  if (length(slide) != 1 || !is.numeric(slide))
    throw("Argument 'slide' must be an integer: ", slide, collapse=", ");
  
  if (!is.null(what) && !is.element(what, getFieldNames(this)))
    throw("Argument 'what' is not refering to a known field: ", what);

  include <- which(getInclude(this, include, exclude, slide=slide));
  
  colWhat <- what;
  if (length(col) == 1) {
    color <- col;
    col <- rep(NA, length.out=nbrOfSpots(this));
    if (!is.numeric(color) && substring(color,1,1) != "#" && !is.element(color, colors())) {
      rg <- getForeground(this, slides=slide);
      R <- rg$R;
      G <- rg$G;
      rm(rg);
      cols <- MicroarrayData$createColors(R,G, type="RG", palette=palette, A.range=A.range, M.range=M.range);
      col[include] <- cols[include];
      rm(cols);
      # Generate the colors automatically using the value of specified field.
      # col[include] <- getColors(this, what=colWhat, slide=slide, include=include, palette=color, log=log);
    } else {
      # Ex) col="blue" etc gets here
      col[include] <- color;
    }
  }

  xy <- getSpotPosition(this, slides=slide);
  # Upper *left* corner is (0,0)
  x <- xy$x; y <- -xy$y; # Will also be used as xlab and ylab.
  
  if (missing(yaxt))
    yaxt0 <- "n";

  Device$setPlotParameters(object=this, fcn="plotSpatial", what=what, slide=slide, look=look);

  if (is.null(cmd)) {
    plot(x,y, pch=pch, col=col, cex=cex, yaxt=yaxt0, ...);
  
    if (missing(yaxt)) {
      # Add the y-axis with the correct ticks
      yaxp <- par("yaxp");
      yaxis <- seq(yaxp[1],yaxp[2], length=yaxp[3]+1);
      axis(2, at=yaxis, labels=-yaxis);
    }
  } else if (cmd == "points") {
    col <- col[include];
    points(x[include],y[include], cex=cex, col=col, pch=pch, ...);
  } else if (cmd == "highlight") {
    col <- col[include];
    points(x[include],y[include], cex=cex, col=col, pch=pch, ...);
  } else if (cmd == "text") {
    col <- col[include];
    text(x[include],y[include], cex=cex, col=col, ...);
  }
})


setMethodS3("plotSpatial3d", "GenePixData", function(this, field="Log Ratio", ...) {
  NextMethod("plotSpatial3d", this, field=field, ...);
})



setMethodS3("normalizeGenewise", "GenePixData", function(this, fields=NULL, bias=0, scale=1, ...) {
  if (is.null(fields)) {
    fields <- getFieldNames(this);
    exclude <- c("Block", "Column", "Row", "Name", "ID", "X", "Y", "Flags");
    fields <- setdiff(fields, exclude);
  }
  NextMethod("normalizeGenewise", this, fields=fields, bias=bias, scale=scale, ...);
})




setMethodS3("getBackground", "GenePixData", function(this, which=c("median","mean"), slides=NULL, channels=1:2, ...) {
  which <- match.arg(which);
  slides <- validateArgumentSlides(this, slides=slides);
  channels <- validateArgumentChannels(this, channels=channels);

  if (which == "mean") {
    fields <- c("B<channel[1]> Mean", "B<channel[2]> Mean")
  } else if (which == "median") {
    fields <- c("B<channel[1]> Median", "B<channel[2]> Median")
  }
  fields <- convertFieldNames(this, fields, channels=channels);

  # Assert that the fields do really exist.
  if (!all(fields %in% getFields(this))) {
    throw("The background estimates ", paste(fields, collapse=" and "), 
          " is not part of this ", data.class(this), " object.");
  }

  RGData(R=this[[fields[1]]][,slides], G=this[[fields[2]]][,slides], layout=getLayout(this));
})



setMethodS3("getForeground", "GenePixData", function(this, which=c("median", "mean"), slides=NULL, channels=1:2, ...) {
  which <- match.arg(which);
  slides <- validateArgumentSlides(this, slides=slides);
  channels <- validateArgumentChannels(this, channels=channels);

  if (which == "mean") {
    fields <- c("F<channel[1]> Mean", "F<channel[2]> Mean")
  } else if (which == "median") {
    fields <- c("F<channel[1]> Median", "F<channel[2]> Median")
  }
  fields <- convertFieldNames(this, fields, channels=channels);

  # Assert that the fields do really exist.
  if (!all(fields %in% getFields(this))) {
    throw("The foreground estimates ", paste(fields, collapse=" and "), 
          " is not part of this ", data.class(this), " object.");
  }

  RGData(R=this[[fields[1]]][,slides], G=this[[fields[2]]][,slides], layout=getLayout(this));
})





setMethodS3("anonymize", "GenePixData", function(this, ...) {
  # Anonymize a copy of the Layout object; others might use this one.
  layout <- clone(getLayout(this));
  anonymize(layout, ...);
  
  # Assume the same layout on all slides!
  slides <- 1:nbrOfSlides(this);
  this[["Name"]][,slides] <- getName(layout);
  this[["ID"]][,slides] <- getID(layout);
})


# Wrapper for nasty field name 'Rgn R\262' which is 'Rgn R^2'.
setMethodS3("getRgn R2", "GenePixData", function(this, ...) {
  as.matrix(this[["Rgn R\262"]]);
})



#########################################################################/**
# @RdocMethod getForegroundSD
#
# @title "Gets the standard deviation of the foreground pixels"
#
# @synopsis
#
# \description{
#  @get "title". 
# }
#
# \value{
#   Returns a list of matrices that contain the estimated standard deviation
#   of the pixels in the foreground region of the spots.
# }
#
# \details{
#   The GenePix fields \code{"F532 SD"} (green) and \code{"F635 SD"} (red) 
#   are returned.
# }
#
# \examples{
#   gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#   raw <- getRawData(gpr)
#   sd <- getForegroundSD(gpr)
#   raw$RSD <- sd$RSD; raw$GSD <- sd$GSD;
#   sd <- getBackgroundSD(gpr)
#   raw$RbSD <- sd$RbSD; raw$GbSD <- sd$GbSD;
#   subplots(4)
#   plot(raw, "RSDvsR", col="red")
#   plot(raw, "RbSDvsRb", col="red")
#   plot(raw, "GSDvsG", col="green")
#   plot(raw, "GbSDvsGb", col="green")
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################  
setMethodS3("getForegroundSD", "GenePixData", function(this, slides=NULL, channels=1:2, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  channels <- validateArgumentChannels(this, channels=channels);

  fields <- c("F<channel[1]> SD", "F<channel[1]> SD");
  fields <- convertFieldNames(this, fields, channels=channels);
  RSD <- this[[fields[1]]][,slides];
  GSD <- this[[fields[2]]][,slides];
  list(RSD=as.matrix(RSD), GSD=as.matrix(GSD));
})



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
setMethodS3("getForegroundSE", "GenePixData", function(this, slides=NULL, channels=1:2, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  channels <- validateArgumentChannels(this, channels=channels);

  sd <- getForegroundSD(this, slides=slides);
  n <- getArea(this, slides=slides);
  sd <- lapply(sd, FUN=function(x) x / n);
  names(sd) <- gsub("SD", "SE", names(sd));

  sd;
})


#########################################################################/**
# @RdocMethod getBackgroundSD
#
# @title "Gets the standard deviation of the background pixels"
#
# @synopsis
#
# \description{
#  @get "title". 
# }
#
# \value{
#   Returns a list of matrices that contain the estimated standard deviation
#   of the pixels in the background region of the spots.
# }
#
# \details{
#   The GenePix fields \code{"B532 SD"} (green) and \code{"B635 SD"} (red) 
#   are returned.
# }
#
# \examples{\dontrun{See help for getForegroundSD() for this class.}}
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################  
setMethodS3("getBackgroundSD", "GenePixData", function(this, slides=NULL, channels=1:2, ...) {
  slides <- validateArgumentSlides(this, slides=slides);
  channels <- validateArgumentChannels(this, channels=channels);

  fields <- c("B<channel[1]> SD", "B<channel[1]> SD");
  fields <- convertFieldNames(this, fields, channels=channels); 
  RbSD <- this[[fields[1]]][,slides];
  GbSD <- this[[fields[2]]][,slides];
  list(RbSD=as.matrix(RbSD), GbSD=as.matrix(GbSD));
})


############################################################################
#
#   S c a n n e r   S e t t i n g s   e t c .
#
############################################################################
setMethodS3("getWavelengths", "GenePixData", function(this, ...) {
  value <- getHeaderField(this, "^Wavelengths.*");
  if (is.null(value)) {
    value <- matrix(c(635, 532), nrow=2, ncol=nbrOfSlides(this));
    return(value);
  }

  # Remove "nm" if that is used.
  value <- gsub("[^0-9]", "", value);
  matrix(as.integer(value), ncol=nbrOfSlides(this));
})


setMethodS3("getPMTGain", "GenePixData", function(this, ...) {
  value <- getHeaderField(this, "^PMT.*");
  if (is.null(value))
    return(NULL);
  matrix(as.integer(value), ncol=nbrOfSlides(this));
})

setMethodS3("getLaserPower", "GenePixData", function(this, ...) {
  value <- getHeaderField(this, "^LaserPower.*");
  if (is.null(value))
    return(NULL);

  matrix(as.double(value), ncol=nbrOfSlides(this));
})

setMethodS3("getLaserOnTime", "GenePixData", function(this, ...) {
  value <- getHeaderField(this, "^LaserOnTime.*");
  if (is.null(value))
    return(NULL);
  matrix(as.double(value), ncol=nbrOfSlides(this));
})

setMethodS3("getLinesAveraged", "GenePixData", function(this, ...) {
  value <- getHeaderField(this, "^LinesAveraged.*");
  if (is.null(value))
    return(1);
  matrix(as.integer(value), ncol=nbrOfSlides(this));
})

setMethodS3("getTemperature", "GenePixData", function(this, ...) {
  value <- getHeaderField(this, "^Temperature.*");
  if (is.null(value))
    return(NULL);
  matrix(as.double(value), ncol=nbrOfSlides(this));
})

setMethodS3("getFocusPosition", "GenePixData", function(this, ...) {
  value <- getHeaderField(this, "^FocusPosition.*");
  if (is.null(value))
    return(NULL);
  matrix(as.double(value), ncol=nbrOfSlides(this));
})

setMethodS3("getPixelSize", "GenePixData", function(this, ...) {
  value <- getHeaderField(this, "^PixelSize.*");
  if (is.null(value))
    return(NULL);
  matrix(as.integer(value), ncol=nbrOfSlides(this));
})


############################################################################
#
#   A r r a y   d i m e n s i o n s
#
############################################################################
setMethodS3("getArrayLeft", "GenePixData", function(this, slide=NULL, ...) {
  slide <- validateArgumentSlide(this, slide=slide);
  getLeftEdge(getSpotPosition(this, slide=slide));
})

setMethodS3("getArrayRight", "GenePixData", function(this, slide=NULL, ...) {
  slide <- validateArgumentSlide(this, slide=slide);
  getRightEdge(getSpotPosition(this, slide=slide));
})

setMethodS3("getArrayTop", "GenePixData", function(this, slide=NULL, ...) {
  slide <- validateArgumentSlide(this, slide=slide);
  getTopEdge(getSpotPosition(this, slide=slide));
})

setMethodS3("getArrayBottom", "GenePixData", function(this, slide=NULL, ...) {
  slide <- validateArgumentSlide(this, slide=slide);
  getBottomEdge(getSpotPosition(this, slide=slide));
})

setMethodS3("getArrayWidth", "GenePixData", function(this, slide=NULL, ...) {
  getMaxWidth(getSpotPosition(this, slide=slide));
})

setMethodS3("getArrayHeight", "GenePixData", function(this, slide=NULL, ...) {
  getMaxHeight(getSpotPosition(this, slide=slide));
})

setMethodS3("getArrayAspectRatio", "GenePixData", function(this, slide=NULL, ...) {
  getAspectRatio(getSpotPosition(this, slide=slide));
})



############################################################################
#
#   S p o t    p r o p e r t i e s  (non-channel specific)
#
############################################################################

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
#   \item{slides}{Specifying for which slides the spot positions should
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
# \examples{
#   gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))
#
#   # Gets the spot positions
#   xy <- getSpotPosition(gpr)
# }
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getSpotPosition", "GenePixData", function(this, slides=NULL, index=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (any(index < 1) || any(index > nbrOfSpots(this)))
    throw("Argument 'index' is out of range.");
  
  if (is.null(index))
    index <- 1:nbrOfSpots(this);

  xField <- "X";
  yField <- "Y";
  if (!hasField(this, xField) || !hasField(this, yField)) {
    throw("This ", data.class(this), " object is missing the fields '", xField, "' and/or '", yField, "', which specifies the physical position of the spots.");
  }

  x <- this[[xField]][index,slides];
  y <- this[[yField]][index,slides];
  SpotPosition(x=x, y=y);
})




#########################################################################/**
# @RdocMethod getArea
#
# @title "Gets the area of the spot"
#
# \description{
#  @get "title". The area is the number of pixels identified (segmentated)
#  to be within the spot limits.
# }
#
# @synopsis
#
# \arguments{
#   \item{slides}{Specifying for which slides the spot area should
#     be extracted. If @NULL, all slides are considered.}
#   \item{include}{The spots for which the spot area is returned.
#       If @NULL all spots are considered.}
# }
#
# \value{
#   Returns a @matrix.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getArea", "GenePixData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  # "'F Pixels' - the total number of feature pixels."
  as.matrix(this[["F Pixels"]][include,slides]);
})



setMethodS3("getBackgroundArea", "GenePixData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  # "'B Pixels' - the total number of background pixels."
  as.matrix(this[["B Pixels"]][include,slides]);
})


setMethodS3("getBgArea", "GenePixData", function(this, ...) {
  getBackgroundArea(this, ...);
}, protected=TRUE)




setMethodS3("getCircularity", "GenePixData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  if (hasField(this, "Circularity")) {
    return( as.matrix(this$Circularity[include, slides]) );
  }

  # "'F Pixels' - the total number of feature pixels."
  area <- this[["F Pixels"]][include,slides];
  # "'Dia.' - the diameter in µm of the feature-indicator."
  dia <- this[["Dia."]][include,slides];
  # "'PixelSize=10' - Resolution of each pixel in µm."
  # (assume same pixel size for all pixels)
  pixelSize <- getPixelSize(this)[,slides];
  # Diameter is pixels.
  d <- dia / pixelSize;
  # area = pi*r^2 = pi*(d/2)^2 = pi/4 * d^2  <=> pi/4 = area / d^2
  # However, this is only true for circles. So calculate
  #   circularity = area / (pi/4 * d^2)
  # where 0 <= circularity <= 1, but for discretization reasons it
  # might be larger than one too.
  circularity <- area / (pi/4 * d^2)
  as.matrix(circularity);
})


setMethodS3("getAbscent", "GenePixData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  # "'Flags' - the type of flag associated with a feature."
  # A value of -100 means a "bad" spot.
  # A value of  -75 means a "abscent" spot.
  # A value of  -50 means a "not found" spot.
  as.matrix(this[["Flags"]][include,slides] == -75);
})

setMethodS3("getBad", "GenePixData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  # "'Flags' - the type of flag associated with a feature."
  # A value of -100 means a "bad" spot.
  # A value of  -75 means a "abscent" spot.
  # A value of  -50 means a "not found" spot.
  as.matrix(this[["Flags"]][include,slides] == -100);
})

setMethodS3("getNotFound", "GenePixData", function(this, slides=NULL, include=NULL, ...) {
  slides <- validateArgumentSlides(this, slides=slides);

  if (is.null(include))
    include <- seq(nbrOfSpots(this));

  # "'Flags' - the type of flag associated with a feature."
  # A value of -100 means a "bad" spot.
  # A value of  -75 means a "abscent" spot.
  # A value of  -50 means a "not found" spot.
  as.matrix(this[["Flags"]][include,slides] == -50);
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#   Channel names etc.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


#########################################################################/**
# @RdocMethod convertFieldNames
#
# @title "Convert field names that contains patterns"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{patterns}{A @vector of @character strings containing field names
#     with patterns to be converted.}
#   \item{channels}{An @integer @vector.}
#   \item{channelNames}{A @character @vector.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns the @vector of @character string of the converted field names.
# }
#
# @author
#
# \examples{
#   print(GenePixData$convertFieldNames(c("area", "F<channel[1]> Mean")));
# }
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################s
setMethodS3("convertFieldNames", "GenePixData", function(this, patterns, channels=NULL, channelNames=NULL, ...) {
  if (length(patterns) == 0)
    return(patterns);

  if (is.null(channelNames))
    channelNames <- getChannelNames(this);

  if (!is.null(channels)) {
    if (!is.numeric(channels))
      throw(InternalErrorException("Argument 'channels' must be numeric or NULL."));
    if (any(channels < 1 || channels > length(channelNames)))
      throw(InternalErrorException("Argument 'channels' is out of range."));
    channelNames <- channelNames[channels];
  }

  for (kk in seq(along=channelNames)) {
    search <- sprintf("<channel\\[%d\\]>", kk);
    patterns <- gsub(search, channelNames[kk], patterns);
  }

  patterns;
}, private=TRUE)



#########################################################################/**
# @RdocMethod getChannelNames
#
# @title "Gets the names of the channels"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{digitsOnly}{If @TRUE, all non-digit characters are excluded from
#     the channel names.}
#   \item{settings}{Internal use only.}
# }
#
# \value{
#  Returns an @character string @vector. An element with value @NA is a
#  channel that had a zero-length name.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################s
setMethodS3("getChannelNames", "GenePixData", function(this, digitsOnly=TRUE, settings=NULL, ...) {
  # 1. Get type of GPR file
  type <- NULL;
  if (!is.null(settings)) {
    idx <- which(names(settings) == "Type");
    if (length(idx) == 1) 
      type <- settings[idx];
  }
  if (is.null(type))
    type <- getHeaderField(this, "^Type");
  if (is.null(type))
    type <- "-1";

  # Get the version number
  type <- gsub("GenePix Results", "", type);
  type <- as.double(type);

  wavelengths <- NULL;

  # 2. Check names of wavelengths in the 'settings' argument....
  if (!is.null(settings)) {
    idx <- which(names(settings) == "Wavelengths");
    if (length(idx) == 1) 
      wavelengths <- settings[idx];
  }

  # 3. If not found there, check in the header fields...
  if (is.null(wavelengths)) {
    # Get the number of possible channels
    wavelengths <- getHeaderField(this, "^Wavelengths.*");
  }

  if (!is.null(wavelengths)) {
    # Split between tabs
    wavelengths <- unlist(strsplit(as.character(wavelengths), split="\t"));
  }



  # 3. If not found there either, go for default names (depending on version)
  nbrOfWavelengths <- length(wavelengths);
  if (nbrOfWavelengths == 0) {
    if (type < 2) {
      wavelengths <- c("635", "532");
    } else {
      wavelengths <- c("635", "532", "3", "4");
    }
  } else if (nbrOfWavelengths < 4) {
    wavelengths <- c(wavelengths, seq(from=4-nbrOfWavelengths+1, to=4));
  }

  # 4. Keep only digits, i.e. exclude all non-digit characters
  #    from the wavelength strings.
  #    Comment: Although the user interface of Axon GenePix Pro only allows
  #    you to enter three-digit labels for each wavelength/channel, for some reason
  #    does Axon's online GPR example files contain examples with
  #      "Wavelengths=635 nm\t532 nm" (v5.0.0.x)
  #      "Wavelengths=635 nm\t532 nm\tn/a\tn/a" (v4.1.1.x)
  if (digitsOnly) {
    wavelengths <- gsub("[^0-9]", "", wavelengths);
  }

  # 5. Set empty channel names (from for instance "n/a" strings etc) to NA
  wavelengths[nchar(wavelengths) == 0] <- NA;

  wavelengths;
}, private=TRUE);


#########################################################################/**
# @RdocMethod validateArgumentChannels
#
# @title "Validates the argument 'channels'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{channels}{Either an @integer or a @character string @vector.}
#   \item{minLength, maxLength}{The minimum and the maximum length of the
#     argument \code{channels}. These arguments are used internally.}
# }
#
# \value{
#  Returns an @integer @vector of the same length as \code{channels} of
#  channel indicies that are known to be valid.
#  If the argument was invalid an @see "R.oo::Exception" is thrown.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("validateArgumentChannels", "GenePixData", function(this, channels, minLength=length(channels), maxLength=minLength, ...) {
  if (is.null(channels))
    throw("Argument 'channels' is NULL. It must be a vector.");
  if (!is.vector(channels)) 
    throw("Argument 'channels' is not a vector.");
  if (length(channels) < minLength)
    throw("Argument 'channels' must contain at least ", minLength, ": ", length(channels));
  if (length(channels) > maxLength)
    throw("Argument 'channels' must not contain more than ", maxLength, ": ", length(channels));

  # Get the known channel names (without non-digit characters).
  channelNames <- getChannelNames(this);
  nchannels <- length(channelNames);
  if (length(nchannels) == 0)
    throw(InternalErrorException("GPR file format error. The 'Wavelengths' header exist but is empty!"));

  # Treat the special cases where 'channels' is in "R", "G", "Cy3" and "Cy5".
  if (is.character(channels)) {
    channels[channels == "R"] <- 1;
    channels[channels == "G"] <- 2;
    channels[channels == "Cy5"] <- 1;
    channels[channels == "Cy3"] <- 2;
  }
  
  # If numeric, consider it as indicies.
  if (is.numeric(channels)) {
    if (any(channels < 1 | channels > nchannels)) {
      throw("Channel index out of range [1,", nchannels, "]: ", paste(channels, collapse=", "));
    }
    return(channels);
  }

  # Now, assume everything else can be treated as characters.
  channels <- as.character(channels);

  # Match the requested channels with the known ones.
  pos <- match(channels, channelNames);
  if (any(is.na(pos))) {
    throw("Some of the specified channels are unknown: ", 
              paste(channels[channels(is.na(pos))], collapse=", "));
  }

  pos;
}, private=TRUE)

############################################################################
# HISTORY:
# 2006-05-04
# o BUG FIX: An error would be thrown in getChannelNames() if the GPR header
#   contained two headers with names containing "Type".  Thanks Valtteri 
#   Wirta, KTH for reporting this.
# 2006-02-23
# o Appended '...' to all methods.
# o readOneFile() now tries to read the results table using different
#   kind of quotes.  This was necessary because one user send me a file
#   where the names contains newlines.
# o BUG FIX: convertFieldNames() would not recognize channel 3 and 4.
# 2005-10-21
# o Replace 'overwrite' arguments with 'mustNotExist' in calls to Arguments. 
# 2005-07-19
# o Replaced all path="" arguments to path=NULL.
# 2005-06-11
# o Making use of Arguments in R.utils.
# 2005-05-03
# o Updated regular expressions.
# 2004-08-17
# o Updated the error message of getForeground(), because it refered to 
#   the "background".
# 2004-08-15
# o Renamed as.RawData() to getRawData().
# 2004-05-06
# o BUG FIX: GenePixData$read() opened a file/connection without closing it
#   eventually resulting in an error with too many opened connections. The
#   was a single left over code of line that opened a file and should not
#   be there.
# 2004-04-26
# o Recoded internally to make use of new methods of class SpotPosition.
# 2004-01-22
# o Added the readHeader() method.
# 2004-01-21
# o No update, but Gordon Smyth reported that using comment.char="" halfs
#   the reading time compared to default comment.char="#". "Unfortunately",
#   this was already done for other reasons. See below.
# 2004-01-14
# o Added argument 'slides' and 'channels' to all methods where applicable.
#   Another big step towards multi-channel data analysis.
# o Added getChannelNames() and validateArgumentChannels().
# o Update getHeaderField() to return NULL if no header fields exist.
# 2004-01-13
# o Now read() and write() of GenePixData is sensitive to the "Wavelengths"
#   header field, if it exists. Its contents is parsed (basically string
#   split around '\t' and removes all known non-digits) and used as 
#   channel labels. This replaces the hard coded "635", "532", "3" and "4"
#   as used before. Another step towards multi-color/channel analysis.
# 2003-12-01
# o BUG FIX: In append() I did rep(NA, nrow=...) instead of 
#   rep(NA, length=...). Thanks to Jarmo Schrader at Umea University, 
#   Sweden, for reporting this.
# 2003-11-18
# o When reading GPR files the estimation/guessing of the number of grid 
#   rows was improved. Before it happend that it failed because some spot
#   positions (x,y) were NAs. If that occurs, it will try again with 
#   another set of spots and so on.
# o When reading GPR files with "missing" spots, i.e. when the printtip
#   groups are not completly rectangular, "missing" spots are inserted at
#   the missing (Block,Row,Column) with all other values set to NA. The
#   speed for doing this has been improved dramatically.
# 2003-10-30
# o Now using getSpotPosition() whereever possible. This will make the code
#   more robust for future changes in field names.
# o Added getArrayXXX() method to get physical dimensions, margins etc.
#   Plans to use this for automatic aspect ratio when calling plotSpatial().
# o Added support for highlight(), points() and text() after doing a
#   plotSpatial() on a GenePixData object.
# 2003-10-26
# o Added support to read *.gpr.gz files too.
# 2003-10-19
# o getForeground() and getBackground() now returns the *median* by default.
#   as.RawData() was already update like a few months ago.
# 2003-09-29
# o BUG FIX: append() did not append header fields correctly after the new 
#   update. Now we assume that all header fields can be stored as column
#   vectors, i.e. one column for each slide. It is not clear whether this
#   is true for all GPR formats or not, but we assume so.
# 2003-09-24
# o Updated append() to be more robust against different GPR version.
# 2003-08-29
# o Created private and static method getKnownHeadersWithType() to be used
#   by both the read- and the write-() methods. The source of that method
#   is also much more clear about what columns are introduced in what
#   GPR file version.
# o Added match.arg() where applicable.
# 2003-07-29
# o Updated as.RawData() to use the *median* estimates instead of the *mean*
#   estimates of the foreground and background signals.
# 2003-06-23
# o Added getBackgroundArea().
# o Added getWavelengths(), getPMTGain(), getLaserPower(), getLaserOnTime(),
#   getLinesAveraged(), getTemperature(), getFocusPosition(), getPixelSize().
# o Made the version field into a matrix too.
# o Made all GPR header fields into regular fields too. This means that now
#   the header fields are stored for *all* slides if several slides are
#   appended together. Moreover, if a header field contains TABS ("\t") 
#   they are used to split the field into a vector of string value to make
#   it easier for the user to retrieve for instance the wavelength for
#   each channel etc. Old 'settings' field is removed.
# o Updated append() to support the header fields.
# o Added protected getHeaderFields(), getHeaderField() & setHeaderField().
# 2003-06-15
# o Added getForegroundSE().
# o Added getForegroundSD() and getBackgroundSD().
# 2003-04-11
# o Added getBackground() and getForeground().
# 2002-12-24
# o Added plotSpatial3d().
# o When reading GPR files and the gene name field contains "#" as part of
#   the gene name, the internal read.table() treats it as a comment 
#   character and skips the rest of the line. Added comment.char="" (to turn
#   off the interpretation of comments altogether). Thanks Sangsoo, South
#   Korea, for the report and the fix.
# 2002-12-21
# o When reading slides from files each slide is now named as the filename.
# 2002-12-05
# o BUG FIX: The get*() methods did not return a matrix if there was only
#   one slide.
# 2002-11-24
# o Now read() can also read data files where some spots are missing. This
#   is done by reconstructing the missing rows and adding NA's.
# 2002-11-20
# o Minor update to write() of the GenePixData class: Header rows, e.g.
#   "type=value" are not prepended by a spaces anymore.
# 2002-10-04
# o BUG FIX: Updated the read() method to be more forgiving and also accept
#   unknown headers by assuming colClass=NA and not expecting any specific
#   order of the headers. 
# 2002-09-24
# o Changed the attribute 'path="."' to 'path=""' in read().
# o Update the Rdoc's for as.RawData().
# 2002-09-21
# o read() can now read one or several files specified by names or by a
#   pattern. This is identical to readAll(), which is now just calling
#   read() for backward compatibilities.
# 2002-08-20
# o Replace 'append(super(this),other)' with NextMethod("append") in
#   method append().
# 2002-08-06
# o BUG FIX: In previous versions of GenePixData I did not have access to
#   a GenePix v4.0 GPR file and I just wrote up the read/write code from
#   the specifications. Now it has been tested on a real data set.
# o Made the internal code of read() a little bit more general when it comes
#   to converting factor levels to the right colClasses type.
# o Update write() to be a much more memory efficient by retrieving and
#   writing chunks of data instead of all data at once. If this idea turns
#   out to be really great the Spot, ImaGene and ScanAlyze classes will 
#   also be updated.
# o Updated read() to deal with data header rows that contains trailing
#   blanks. Currently I do not know how to deal with leading blanks.
# 2002-07-02
# o BUG FIX: The field 'Log Ratio' was incorrectly translated to levels
#   instead of doubles.
# o BUG FIX: read() failed to read column names that contains
#   non-ASCII 0-127 characters, e.g. "Rgn R^2" where '^2' is ASCII
#   hex 0xB2, octal \262, or decimal 178.
#   Added a "getRgn R2"() wrapper for the "Rgn R^2" field, which allows
#   the field to be access as gpr[["Rgn R2"]] :)
# o Updated write() to work for more different kinds of versions.
# o Updated read() internally so it is easier to update to recognize new
#   GPR file header versions. Added code for GPR file version v3.0, but
#   haven't been able to test it.
# 2002-07-01
# o BUG FIX: NextMethod() should only take arguments method name and object
#   and *nothing* more. It generates an error.
# o BUG FIX: Memory leak in anonymize(). Forgot to delete temporary Layout
#   object.
# 2002-06-24
# o BUG FIX: The read()/write() fix of 2002-05-20 (see below) was not fully
#   correct for write(), which remained broken. A small update did it!
# 2002-05-20
# o BUG FIX: When started to use [R] v1.5.0 read() and write() would give a
#   cohersion error: No method or default for coercing "character" to "NA".
#   This was due to use "NA"'s instead of NA's in colClasses specifications
#   to read.table().
# 2002-05-03
# o BUG FIX: GenePixData$read() could not read Demo.gpr from Axon. It was
#   a too old version of file format. Can read now!
# o Added anonymize().
# o Added getBad(), getAbscent(), getNotFound().
# o Added Rdoc comments about the Flags field; -100, -75, -50 & 0 values.
# o Added getArea(), getCircularity() and getBgArea().
# 2002-04-05
# * Renamed getSpotPositions() to getSpotPosition() and added the argument
#   index=NULL to it.
# 2002-04-03
# * Added argument 'solve=FALSE' to read().
# * read() and write() now supports File objects.
# 2002-03-29
# * BUG FIX: Added as.character() to getID() and getName() just to make sure
#   it is not a factor that is returned.
# 2002-03-25
# * Added getSpotPositions() and plotSpatial().
# 2002-02-28
# * Made class implements interface Morpable.
# 2002-02-27
# * Added strip.white=TRUE too all scan's.
# 2002-02-26
# * Removed the dye.swap argument in as.RawData().
# * Added Rdoc comments for write. Updated the other Rdoc comments.
# * write() now adds (if non-exists) a Creator field.
# 2002-02-25
# * Can now write GenePixResult data format files.
# * Correction in reading the 'settings' from a gpr file; quotes are now
#   excluded.
# * GenePixData now stores each field as an object field of type matrix.
#   For this reason the class ResultsData becomes obsolete.
# * Updated with setMethodS3().
# 2002-01-19
# * In read() colClasses are now set to increase the speed of read.table().
# * "BUG FIX": Started to use regexpr() to compare the read header with
#   known headers. This takes care of the annoying warning about "Rgn R2"
#   on Windows systems.
# 2001-11-17
# * Updated readAll() to also include pattern matching.
# 2001-11-12
# * Bug fix: Used version instead of self$version in read.GenePixData().
# 2001-07-12
# * Bug fix: Forgot to do quote="" in all scan() and read.table() calls.
#   Trying to read GPR files with cells including "'" would give an error.
# 2001-07-11
# * Removed a lot of the fields that weren't much of interest for public 
#   use anyway.
# * Removed the updateLayout method and put everything into read.
# * Renamed the class from GenePixResultsData to SpotData.
# 2001-07-06
# * Removed getSpotAnnotations().
# 2001-06-30
# * Now unfactorizes "Name" and "ID" and "Log.Ratio" in read().
# 2001-06-29
# * Added updateLayout() and removed getLayout(), which is now in the
#   MicroarrayData class.
# 2001-06-09
# * Improved the Rdoc comments.
# 2001-06-07
# * Forgot to make read() static which gave strange behavior, especially
#   since I was showing everything down at Ernest Gallo Research Lab.
# 2001-04-01
# * Added method getSpotAnnotations().
# 2001-03-27
# * Removed getRGData(); use getRawData()$getSignal() or 
#   getRawData()$getSignal(bg.subtract=TRUE) instead.
# 2001-03-11
# * Created from older com.braju.genepix.R. Now made into a class!
############################################################################
