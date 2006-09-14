###########################################################################/**
# @RdocDefault readBaseFile
#
# @title "Low-level function to read a BASE file (stream) from a connection or a file"
#
# \description{
#  @get "title". 
# }
#
# @synopsis
#
# \arguments{
#   \item{con}{A @connection or a @character string filename.}
#   \item{...}{Arguments such as the important \code{extractSpotsData},
#     passed to @see "readBaseFileSection".}
#   \item{verbose}{Either a @logical, a @numeric, or a @see "R.utils::Verbose"
#     object specifying how much verbose/debug information is written to
#     standard output. If a Verbose object, how detailed the information is
#     is specified by the threshold level of the object. If a numeric, the
#     value is used to set the threshold of a new Verbose object. If @TRUE, 
#     the threshold is set to -1 (minimal). If @FALSE, no output is written
#     (and neither is the \link[R.utils:R.utils-package]{R.utils} package required).
#   }
# }
#
# \value{
#   Returns a named @list structure containing the parsed BASE structure.
#   The names of the elements in the list are the same as the names (types)
#   of the sections as given by section header 'section'.
# }
#
# @examples "../incl/readBaseFile.Rex"
#
# \seealso{
#   @see "writeBaseFile".
#   @see "readBaseFileSection".
#   See \code{stdin()} in \link[base]{showConnections} to read from standard input.
# }
#
# @author
#
# \references{
#   [1] BASE - BioArray Software Environment,\cr
#       \url{http://base.thep.lu.se/}\cr
#
#   [2] Carl Troein, How to write a plugin for BASE, April 2003.\cr
#       \url{http://base.thep.lu.se/documentation/development/plugins.txt}\cr
#       \url{http://opensource.microarray.omrf.org/wiki/bin/view/BASE/PluginWritingHowto}
# }
#
# @keyword file
# @keyword IO
#*/###########################################################################
setMethodS3("readBaseFile", "default", function(con, ..., verbose=FALSE) {
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

  readMagic <- function(con, ...) {
    # Should read "BASEfile" on the first line. (Skip empty lines)
    while (TRUE) {
      line <- readBaseFileLines(con, n=1);
      if (length(line) == 0)
        break;
      line <- trim(line);
      if (nchar(line) > 0)
        break;
    }
    
    if (!identical(line, "BASEfile")) {
      lines <- readLines(con, n=10);
      throw("File format error: Not a BASE file. Missing 'BASEfile' on first line: '", line, "'. Next ten (raw) lines are: '", paste(lines, collapse="\n"), "'");
    }

    verbose && cat(verbose, "Identified the required 'BASEfile' header.");

    TRUE;
  } # readMagic()

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

  # Argument 'verbose':
  if (inherits(verbose, "Verbose")) {
  } else if (is.numeric(verbose)) {
    verbose <- Verbose(threshold=verbose);
  } else {
    verbose <- as.logical(verbose);
    if (verbose)
      verbose <- Verbose(threshold=-1);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read BASE file stream
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Try to read the magic first line.
  readMagic(con);

  # Try to read sections
  sections <- list();
  while (TRUE) {
    section <- readBaseFileSection(con, ..., verbose=verbose);
    if (is.null(section))
      break;
    sections <- append(sections, section);
  } # while()

  sections;
})



############################################################################
# HISTORY: 
# 2005-12-20
# o Added more help on what the names of the returned elements are.
# 2005-06-11
# o Now argument 'verbose' can be logical, numeric or a Verbose object.
# 2005-05-31
# o Cleaned up code. Most arguments are now in readBaseFileSection().
# 2005-05-25
# o Moved readBaseFileSection() to its own file.
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
