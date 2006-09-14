###########################################################################/**
# @RdocDefault readBaseFileSection
#
# @title "Low-level function to read a BASE file section from a connection or a file"
#
# \description{
#  @get "title". This a supportive function to readBaseFile().
# }
#
# @synopsis
#
# \arguments{
#   \item{con}{A @connection or a @character string filename.}
#   \item{suppressSingleAssayLabels}{If @TRUE and only a single assay is
#     read, field names in data table is non suffixed with assay name.}
#   \item{extractSpotsData}{If @TRUE, data tables in 'spots' sections are
#     written to seperate files such that each assay has its own file.
#     Files are written as tab-delimited files to current directory, 
#     using filename format "spots\_<assay>.txt".
#     Written data is removed from the final list structure and replaced
#     with a header \code{dataFiles} of the filenames for each assay.}
#   \item{verbose}{Either a @logical, a @numeric, or a @see "R.utils::Verbose"
#     object specifying how much verbose/debug information is written to
#     standard output. If a Verbose object, how detailed the information is
#     is specified by the threshold level of the object. If a numeric, the
#     value is used to set the threshold of a new Verbose object. If @TRUE, 
#     the threshold is set to -1 (minimal). If @FALSE, no output is written.
#   }
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @list structure containing the parsed BASE structure.
# }
#
# \section{Out of memory?}{
#   Some BASE files are very large due to large amount of assay data in
#   'spots' sections.  By setting argument \code{extractSpotsData=TRUE},
#   such data is written to separate files for each assay immediately after
#   being read and thereafter removed from memory.  This means that it is
#   in practice possible to read BASE files with almost any number of 
#   assays, also if they are "passed" via stdin.
# }
#
# \section{Slow performance?}{
#   The number of lines in the data table of a BASE file data section is 
#   either not known on before hand or specified by the \code{count} header.
#   If \code{count} is available, the data table is read approximately 
#   eight times faster than if it is not.
#   To the best of our knowledge, starting with BASE v1.2.17 (released
#   November 2005), header \code{count} is included by default and in 
#   previous versions it is not available at all.
#   If \code{count} is not specified, the table is first read line by line
#   to a temporary file and then re-read by @see "base::read.table". (The 
#   reason for this is that \code{read.table()} otherwise will read to end 
#   of file, and not to the first empty line. Unfortunately, this means 
#   that the reading of data is a bit slow.)
# }
#
# \seealso{
#   @see "readBaseFile".
# }
#
# @author
#
# @keyword file
# @keyword IO
#*/###########################################################################
setMethodS3("readBaseFileSection", "default", function(con, suppressSingleAssayLabels=TRUE, extractSpotsData=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local function definitions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readBaseFileLines <- function(con, ...) {
    lines <- readLines(con, ...);
  
    # Remove optional prefix "#". 
    # All header lines (and the "BASEfile" line and the "%") may start 
    # with the character "#" which will not count as part of the header
    # name. The reason for this behavior is that it makes it possible to
    # feed BASEfiles directly into e.g. gnuplot.
    lines <- gsub("^#", "", lines);
  
    # Strip trailing whitespaces
    lines <- gsub("[ \t]*$", "", lines);
  
    lines;
  }

  readBaseFileHeaders <- function(con, ...) {
    readHeader <- function(con, skipEmptyLines=FALSE, ...) {
      verbose && enter(verbose, level=-10, "Reading header");
      on.exit(verbose && exit(verbose), append=TRUE);
  
      # Each header line consists of a name, optionally followed by a
      # tab and a value (which may contain tabs).
      while (TRUE) {
        line <- readBaseFileLines(con, n=1);
        if (length(line) == 0)
          return(NULL);
  
        if (nchar(line) == 0) {
          # Sections are seperated by one or several empty lines.
          # An empty line here => end of section and therefore end 
          # of header.
          if (skipEmptyLines) {
            next;
          } else {
            return(NULL);
          }
        }
  
        if (identical(line, "%")) {
          # If a section has data, it's separated from the headers by a 
          # line containing only the character "%". 
          # A '%' line => end of headers and there end of header.
          return(NULL);
        }
          
        line <- strsplit(line, split="\t")[[1]];
        line <- line[nchar(line) > 0];
        name <- line[1];
        if (!is.na(name)) {
          value <- line[-1];
          if (length(value) > 0) {
            # In case field value is "\\t" (!= '\t');
            value <- unlist(strsplit(value, split="\\t", fixed=TRUE));
          }
          if (length(value) == 0)
            value <- "";
          break;
        }
      } # while()
      header <- list(value);
      names(header) <- name;
      header;
    } # readHeader()
  
    verbose && enter(verbose, "Reading headers");
    on.exit(verbose && exit(verbose), append=TRUE);

    # Read headers
    headers <- list();
    while(TRUE) {
      header <- readHeader(con, skipEmptyLines=(length(headers) == 0));
      if (is.null(header))
        break;
      headers <- c(headers, header);
    }

    headers;
  } # readBaseFileHeaders()
  
  
  
  readBaseFileData <- function(con, nrows=-1, assayFieldSep="_to_", ...) {
    verbose && enter(verbose, "Reading data");
    on.exit(verbose && exit(verbose), append=TRUE);
  
    # Make use of a temporary file to read the data part.
    # Copy all lines until an empty line or EOF is detected.
    tmpfile <- tempfile();
    verbose && cat(verbose, "Copying data lines to temporary file: ", tmpfile);
    fh <- file(tmpfile, open="w");
    on.exit({
      if (!is.null(fh)) close(fh);
      if (file.exists(tmpfile))
        file.remove(tmpfile);
    }, append=TRUE);

    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Alternative 1   
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if (nrows == -1) {
      count <- 0;
      while (TRUE) {
        line <- readLines(con, n=1);
        if (length(line) == 0)
          break;
        line <- gsub("^#", "", line);
  
        # Check for end of data table
        if (nchar(line) == 0)
          break;
  
        # Workaround: If a data cell contains the two characters "\t" (not a
        # TAB), read.table(..., sep="\t") will interprete this as a TAB 
        # indeed. Because of this, we replace all TABs with BACKSPACE in
        # the temporary file.
        # read.table() should provide argument 'allowEscapes' and pass it
        # to scan(). /HB 2005-05-24
        line <- gsub("\t", "\b", line);
    
        writeLines(line, con=fh);
        count <- count + 1;
      }
      rm(line);
    } else {
      lines <- readLines(con, n=nrows);
      lines <- gsub("\t", "\b", lines);
      writeLines(lines, con=fh);
      rm(lines);
      count <- nrows;
    }
  
    verbose && cat(verbose, "Copied ", count, " lines.");
  
    close(fh); fh <- NULL;

    if (count == 0) {  
      verbose && cat(verbose, "No data to read.");
      df <- NULL;
    } else {
      verbose && cat(verbose, "Reading copied data.");
    
      df <- read.table(tmpfile, header=FALSE, sep="\b", quote="", as.is=TRUE, comment.char="", na.strings=c("", "NA", "N/A", "Error"), fill=TRUE);
    
      # Convert all "\t" strings to "\\t". See above.
      for (kk in seq(length=ncol(df))) {
        values <- df[,kk];
        if (is.character(values)) {
          values <- gsub("\t", "\\t", values, extended=FALSE, fixed=TRUE);
          df[,kk] <- values;
        }
      }
      
      verbose && cat(verbose, "Read data table: ", nrow(df), "x", ncol(df));
      verbose && str(verbose, df);
    }
  
    df;
  } # readBaseFileData()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'con':
  if (inherits(con, "connection")) {
  } else {
    filename <- as.character(con);

    if (!file.exists(filename))
      throw("File not found: ", filename);

    # Open a file handler
    con <- file(filename, open="r");
    on.exit(close(con), append=TRUE);
  }

  # Argument 'extractSpotsData':
  extractSpotsData <- as.logical(extractSpotsData);
  if (extractSpotsData) {
    hooks <- getHook("onExit(readBaseFileSection)");
    hook <- extractBaseFileSpotsData;
    setHook("onExit(readBaseFileSection)", hook, action="append");
    on.exit({
      setHook("onExit(readBaseFileSection)", hooks, action="replace");
    }, add=TRUE);
  }

  # Argument 'verbose':
  if (inherits(verbose, "Verbose")) {
  } else if (is.numeric(verbose)) {
    require(R.utils) || throw("Package not available: R.utils");
    verbose <- Verbose(threshold=verbose);
  } else {
    verbose <- as.logical(verbose);
    if (verbose) {
      require(R.utils) || throw("Package not available: R.utils");
      verbose <- Verbose(threshold=-1);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main code
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  assayFieldSep <- "_of_";

  verbose && enter(verbose, "Reading section");
  on.exit(verbose && exit(verbose), append=TRUE);

  # Each section contains a set of headers and an optional set of
  # data lines, which usually have tab-separated columns.

  headers <- readBaseFileHeaders(con);
  if (length(headers) == 0)
    return(NULL);

  sectionName <- headers$section;
  if (is.null(sectionName))
    stop("File format error: Mandatory header 'section' is missing or empty.");

  verbose && cat(verbose, "Section name: ", sectionName);

  # Read optional data
  # If a section has data, it's separated from the headers by a line
  # containing only the character "%".

  if (is.null(headers$count)) {
    nrows <- -1;
  } else {
    nrows <- as.integer(headers$count);
  }
  data <- readBaseFileData(con, nrows=nrows, assayFieldSep=assayFieldSep);

  if (!is.null(data)) {
    # Further processing of BASE structure?
    verbose && enter(verbose, "Trying to add headers to data with ", ncol(data), " unnamed columns");
    columns <- headers$columns;
    if (!is.null(columns)) {
      verbose && cat(verbose, "Identified ", length(columns), " column names: ", paste(columns, collapse=", "));
  
      posAssayData <- which("assayData" == columns);
      if (length(posAssayData) == 1) {
        assayFields <- headers$assayFields;
        nbrOfAssayFields <- length(assayFields);
        assays <- headers$assays;
        nbrOfAssays <- length(assays);
        verbose && enter(verbose, "Creating columns for assay fields");
        verbose && cat(verbose, "Expanding ", nbrOfAssayFields, " assayFields for ", nbrOfAssays, " assay(s)");
  
        if (nbrOfAssays == 1 && suppressSingleAssayLabels) {
          assayDataColumns <- assayFields;
        } else {
          assayDataColumns <- paste(rep(assayFields, times=nbrOfAssays), assayFieldSep, rep(assays, each=nbrOfAssayFields), sep="");
        }
        head <- tail <- NULL;
        if (posAssayData > 1)
          head <- columns[1:(posAssayData-1)];
        if (posAssayData < length(columns))
          tail <- columns[(posAssayData+1):length(columns)];
  
        columns <- c(head, assayDataColumns, tail);
  
        verbose && cat(verbose, "New column names: ", paste(columns, collapse=", "));
        verbose && exit(verbose);
      }
  
      if (ncol(data) > length(columns)) {
        verbose && cat(verbose, "There are more columns in the read data than what is specifed in the header: ", ncol(data), " > ", length(columns));
  
        # Removes (from the right) columns that contains only NA's.
        while (TRUE) {
          values <- data[,ncol(data)];
          if (!all(is.na(values)))
            break;
          data <- data[,-ncol(data)];
          if (ncol(data) <= length(columns))
            break;
        }
      }
  
      if (length(columns) == ncol(data)) {
        verbose && cat(verbose, "Setting column names: ", paste(columns, collapse=", "));
        colnames(data) <- columns;
        columns <- NULL;
      } else { 
        verbose && cat(verbose, "Mismatching number of columns of read data and number of column names: ", ncol(data), " != ", length(columns));
      }
    }
    verbose && exit(verbose);
  } # if (!is.null(data))

  section <- list(headers=headers, data=data);

  res <- callHooks("onExit(readBaseFileSection)", section=section);
  results <- lapply(res, FUN=function(x) x$result);
  for (result in results) {
    if (is.list(result) && "section" %in% names(result)) {
      section <- result$section;
    }
  }

  section <- list(section);
  names(section) <- sectionName;

  section;
}) # readBaseFileSection()



############################################################################
# HISTORY: 
# 2005-12-12
# o Updated the Rdoc comments to mention the new 'counts' header field.
# 2005-07-24
# o Now readBaseFileData() looks for the header 'count', which tells how
#   many rows the following data contains.  This makes it possible to read
#   all of the data at once and not row by row to detect the end of it. 
#   The initial speed up was approximately 8 times!
# 2005-06-26
# o Now the method restores the list of hook function when exiting.
# 2005-06-15
# o BUG FIX: Section with empty data tables gave an error.
# o Added hook 'onExit(readBaseFileSection)'.
# o Added function extractBaseFileSpotsData().
# 2005-06-03
# o Verified for the first that the code can read BASE files with multiple
#   assays in non-serial mode. Serial mode is not implemented yet.
# 2005-05-31
# o Method now strips NA columns in data tables if number of column names
#   is less that the number of read columns.
# 2005-05-25
# o Made readBaseFileSection() public too.
# o Making use of new Verbose class in R.utils.
# 2005-05-24
# o Special treatment of strings "\\t" in data tables is needed! For some
#   reason is read.table() interpreting these as TABs, which is a problem
#   because data tables are in tab-delimited formats. In order to read such
#   files, we first convert all TABs to BACKSPACE and use sep="\b" to read
#   them. In the read data frame, all occurances of '\t' are replaced by
#   "\\t" in each character column.
# o Re-created. Seems to work. The only thing I'm not happy about is that
#   any data parts is read by first copying it line by line into a
#   temporary file in order to detect the end of it (an empty line). 
# 2005-04-10
# o Just a non-working sketch.
############################################################################
