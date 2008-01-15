#########################################################################/**
# @RdocClass ExperimentalSetup
#
# @title "The ExperimentalSetup class"
#
# \description{
#  @classhierarchy
#
#  Represents the experimental setup of cDNA microarray experiments. 
#  An ExperimentalSetup object contains information about such things as
#  samples used, dyes used for each sample on each hybridization, printing
#  information, scanning information, which the image and data files are
#  etc.
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
#
# @author
#
# \examples{
#   \dontrun{
#   # Due to a bug in R CMD check (R v1.7.1) the MicroarrayData$read() call
#   # below will call getKnownSubclasses(), which will generate
#   #   "Error in exists(objectName, mode = "function") : 
#   #	   [2003-07-07 23:32:41] Exception: F used instead of FALSE"
#   # Note that the example still work, just not in R CMD check
#
#   setup <- ExperimentalSetup$read("ScanAlyze.setup", path=system.file("misc", package="aroma"))
#   print(setup)
#   data <- MicroarrayData$read(getDataFiles(setup), path=system.file("misc", package="aroma"))
#   raw <- getRawData(data)
#   ma <- getSignal(raw)
#   subplots(2*nbrOfSlides(ma))
#   for (k in seq(ma)) {
#     plot(ma, slide=k)
#     plotSpatial(ma, slide=k)
#   }
#   }
# }
#*/######################################################################### 
setConstructorS3("ExperimentalSetup", function(...) {
  extend(Object(), "ExperimentalSetup",
    ...
  )
})

setMethodS3("as.character", "ExperimentalSetup", function(this) {
  s <- paste(data.class(this), ": ", nbrOfHybridizations(this), " hybridizations/data files.", sep="");
  s <- paste(s, " Samples used: ", paste(getSampleUsed(this), collapse=", "), ".", sep="");
  s <- paste(s, " Dyes used: ", paste(getDyesUsed(this), collapse=", "), ".", sep="");
  s <- paste(s, " Scanners used: ", paste(getScannersUsed(this), collapse=", "), ".", sep="");
  s <- paste(s, " Image analysis softwares used: ", paste(getImageSoftwaresUsed(this), collapse=", "), ".", sep="");
  s;
})

setMethodS3("print", "ExperimentalSetup", function(x, ...) {
  # To please R CMD check...
  this <- x;

  cat(as.character(this), "\n\n", sep="");
  print(as.data.frame(this, ...))
})

setMethodS3("getDataFiles", "ExperimentalSetup", function(this, channel=NULL) {
  if (is.null(channel))
    channel <- 1:4;
  fields <- paste("raw.file", channel, sep="");
  fields <- intersect(fields, this$.fields);
  res <- list();
  for (field in fields) {
    value <- this[[field]];
    if (!all(is.na(value) | (nchar(value) == 0)))
      res[[field]] <- value;
  }
  res <- as.data.frame(res);
  rownames(res) <- this[["name"]];
  res <- as.matrix(res);
  res;
})

setMethodS3("getLayoutFiles", "ExperimentalSetup", function(this) {
  field <- "print.layout";
  value <- this[[field]];
  if (is.null(value))
    return(NULL);
  res <- list();
  if (!all(is.na(value) | (nchar(value) == 0)))
    res[[field]] <- value;
  res <- as.data.frame(res);
  rownames(res) <- this[["name"]];
  res <- as.matrix(res);
  res;
})

setMethodS3("getSampleUsed", "ExperimentalSetup", function(this) {
  fields <- paste("sample", 1:4, sep="");
  res <- unlist(sapply(fields, FUN=function(field) this[[field]]));
  unique(res[!is.na(res)]);
})

setMethodS3("getDyesUsed", "ExperimentalSetup", function(this) {
  fields <- paste("dye", 1:4, sep="");
  res <- unlist(sapply(fields, FUN=function(field) this[[field]]));
  unique(res[!is.na(res)]);
})

setMethodS3("getScannersUsed", "ExperimentalSetup", function(this) {
  res <- this[["scanner"]];
  unique(res[!is.na(res)]);
})

setMethodS3("getImageSoftwaresUsed", "ExperimentalSetup", function(this) {
  res <- this[["image.software"]];
  unique(res[!is.na(res)]);
})

setMethodS3("nbrOfHybridizations", "ExperimentalSetup", function(this) {
  firstField <- this$.fields[1]
  length(this[[firstField]])
})

setMethodS3("as.data.frame", "ExperimentalSetup", function(x, ...) {
  # To please R CMD check...
  this <- x;

  df <- list()
  for (name in this$.fields)
    df[[name]] <- this[[name]];
  as.data.frame(df, ...);
})


#########################################################################/**
# @RdocMethod read
#
# @title "Reads a file specifying the setup of a microarray experiment"
#
# \description{
#  Static method that read a experiment setup file describing a cDNA 
#  microarray experiment.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the setup file to be read.}
#   \item{path}{Optional path to the setup file.}
#   \item{verbose}{If \@TRUE, information will printed out during
#                  the file reading/parsing.}
# }
#
# \value{
#   Returns an @see "ExperimentalSetup" object.
# }
#
# @author
#
# \examples{
#    setup <- ExperimentalSetup$read("ScanAlyze.setup", path=system.file("misc", package="aroma"))
#    print(setup)
# }
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("read", "ExperimentalSetup", function(static, filename, path=NULL) {
  filename <- Arguments$getReadablePathname(filename, path); 
 
  # Read the setup file
  lines <- readLines(filename);

  # Extract lines that are not blank and are not comment lines.
  incl <- (regexpr("^ *#", lines) == -1);
  incl <- incl & (regexpr("^[ \t]*$", lines) == -1);
  lines <- paste(lines[incl], collapse="\n");

  # Reread the extracted lines using read.table()
  con <- textConnection(lines);
  on.exit(close(con));
  df <- read.table(con, sep="\t", row.names=NULL, header=TRUE, as.is=TRUE, comment.char="", na.strings=c("NA", "-", "?"));

  # Create a ExperimentalSetup object.
  setup <- ExperimentalSetup();
  setup$.fields <- colnames(df);
  for (name in colnames(df))
    setup[[name]] <- df[[name]];

  setup;
}, static=TRUE)



#########################################################################/**
# @RdocMethod findSetupFiles
#
# @title "Locates cDNA microarray experiment setup files"
#
# \description{
#  @get "title" and returns a list of
#  File object, which can be used by the static read() method in this class.
#  By default it searches the current directory, but it can search any number
#  of directories and recursively to an optimal depth.
#
#  Any file with the extension \code{".setup"} is matched. The files are
#  not tried to be read.
# }
#
# @synopsis
#
# \arguments{
#   \item{paths}{Vector of character strings or list of File objects 
#     specifying which directories to be searched.}
#   \item{recursive}{If \@TRUE, all subdirectories will also be search.
#     If \@FALSE, only the specified directory will be searched.
#     If an integer is given, it specifies the maximum depth of the 
#     directory structure to be search; \code{recursive=0} and
#     \code{recursive=FALSE} give the same results.}
# }
#
# \value{
#   Returns a @list of filenames pointing to the setup files found.
# }
#
# @author
#
# \examples{
#    files <- ExperimentalSetup$findSetupFiles(system.file("misc", package="aroma"))
#    print(files)
#    for (file in files) {
#      setup <- ExperimentalSetup$read(file)
#      print(setup)
#    }
# }
#
# \seealso{
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("findSetupFiles", "ExperimentalSetup", function(static, paths, recursive=TRUE) {
  require(R.io) || throw("Package R.io is missing!");

  if (is.logical(recursive) && recursive)
    recursive <- Inf;

  filter <- RegExprFileFilter(".*\\.setup$");
  paths <- lapply(paths, FUN=File);
  found <- list();
  for (path in paths) {
    files <- listFiles(path, filter=filter);
    found <- c(found, files);
    if (recursive > 0) {
      files <- listFiles(path);
      dirs <- files[unlist(lapply(files, FUN=isDirectory))];
      found <- c(found, findSetupFiles(static, dirs, recursive=recursive-1));
    }
  }

  found <- lapply(found, FUN=as.character);

  found;
}, static=TRUE)

############################################################################
# HISTORY:
# 2005-07-19
# o Now returns filenames as character strings.
# 2005-06-11
# o Making use of Arguments in R.utils.
# 2003-05-03
# o Updated all regular expressions.
# 2003-04-26
# o BUG FIX: findSetupFiles() was not declared static.
# 2003-01-07
# o Created!
############################################################################

