###########################################################################/**
# @RdocClass GenericDataFile
#
# @title "The abstract GenericDataFile class"
#
# \description{
#  @classhierarchy
#
#  A GenericDataFile is an object refering to a data file on a file system.
#  Note that this class is abstract and can not be instanciated, but
#  instead you have to use one of the subclasses or the generic 
#  @seemethod "fromFile" method.
# }
# 
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file.}
#   \item{path}{An optional path to the file.}
#   \item{mustExist}{If @TRUE, an exception is thrown if the file does
#     not exists, otherwise not.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# \section{Filename convention}{
#   The filename of an \code{GenericDataFile} is structured as follows:
#   \itemize{
#    \item{filename}{\code{"sample001,a,b,c.CEL"} 
#       (this follows the \R convention (but not the Unix convention)}
#    \item{fullname}{\code{"sample001,a,b,c"}}
#    \item{name}{\code{"sample001"}}
#    \item{tags}{\code{c("a", "b", "c")}}
#    \item{extension}{\code{"CEL"}}
#   }
# }
#
# @author
#
# \seealso{
#   An object of this class is typically part of an @see "GenericDataFileSet".
# }
#*/###########################################################################
setConstructorS3("GenericDataFile", function(filename=NULL, path=NULL, mustExist=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(filename)) {
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=mustExist);
  } else {
    pathname <- NULL;
  }

  # Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  }

  extend(Object(), "GenericDataFile",
    .alias = NULL,
    .pathname = pathname,
    .attributes = list()
  )
}, abstract=TRUE)



setMethodS3("getLabel", "GenericDataFile", function(this, ...) {
  label <- this$label;
  if (is.null(label))
    label <- getName(this, ...);
  label;
}, private=TRUE)

setMethodS3("setLabel", "GenericDataFile", function(this, label, ...) {
  this$label <- label;
  invisible(this);
}, private=TRUE) 


setMethodS3("clone", "GenericDataFile", function(this, clear=TRUE, ...) {
  object <- NextMethod("clone", this, ...);
  if (clear)
    clearCache(object);
  object;
}, private=TRUE)


setMethodS3("equals", "GenericDataFile", function(this, other, ...) {
  if (getPathname(this) == getPathname(other))
    return(TRUE);
  FALSE;
})



###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("as.character", "GenericDataFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  tags <- getTags(this, collapse=",");
  if (!is.null(tags)) {
    s <- c(s, sprintf("Tags: %s", tags));
  }
  s <- c(s, sprintf("Pathname: %s", getPathname(this)));
  s <- c(s, sprintf("File size: %s", getFileSize(this, "units")));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)




###########################################################################/**
# @RdocMethod getPathname
#
# @title "Gets the pathname of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPathname", "GenericDataFile", function(this, ...) {
  this$.pathname;
})




###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path (directory) of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPath", "GenericDataFile", function(this, ...) {
  dirname(this$.pathname);
})




###########################################################################/**
# @RdocMethod getFilename
#
# @title "Gets the filename of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character.
# }
#
# \details{
#  The filename of a file is the pathname excluding any path.
#  For instance, the filename of \code{path/to/foo,a.2,b.ext} is 
#  \code{foo,a.2,b.ext}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFilename", "GenericDataFile", function(this, ...) {
  basename(this$.pathname);
})



###########################################################################/**
# @RdocMethod getFullName
#
# @title "Gets the full name of the file"
#
# \description{
#   @get "title", that is the filename without the extension.
# }
#
# @synopsis
#
# \arguments{
#  \item{aliased}{If @TRUE, and an alias has been set, the alias is 
#     returned, otherwise the default full name is returned.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character.
# }
#
# \details{
#  The full name of a file is the filename excluding any
#  extension (and period).
#  For instance, the full name of \code{path/to/foo,a.2,b.ext} is 
#  \code{foo,a.2,b}.
# }
#
# @author
#
# \seealso{
#   @seemethod "getName".
#   @seemethod "getTags".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFullName", "GenericDataFile", function(this, aliased=FALSE, ...) {
  if (aliased) {
    alias <- getAlias(this);
    if (!is.null(alias))
      return(alias);
  }

  pathname <- this$.pathname;
  if (is.null(pathname))
    return("");

  name <- basename(pathname);

  # Exclude filename extension
  name <- gsub("[.][abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0-9]+$", "", name);

  name;
})


###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character.
# }
#
# \details{
#  The name of a file is the part of the filename without the extension and
#  that preceeds any comma.
#  For instance, the name of \code{path/to/foo,a.2,b.ext} is \code{foo}.
# }
#
# @author
#
# \seealso{
#   @seemethod "getTags".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "GenericDataFile", function(this, ...) {
  name <- getFullName(this, ...);

  # Keep anything before the first comma
  name <- gsub("[,].*$", "", name);
  
  name;
})


setMethodS3("getAlias", "GenericDataFile", function(this, ...) {
  this$.alias;
})

setMethodS3("setAlias", "GenericDataFile", function(this, alias=NULL, ...) {
  if (!is.null(alias)) {
    alias <- Arguments$getFilename(alias);
  }
  
  this$.alias <- alias;
  invisible(this);
})


###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{pattern}{An optional regular expression used to filter out tags.
#     If @NULL, all tags are returned.}
#  \item{collapse}{A @character string used to concatenate the tags. 
#     If @NULL, the tags are not concatenated.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character @vector or @NULL.
# }
#
# \details{
#  The \emph{tags} of a filename are the comma separated parts of the
#  filename that follows the the first comma, if any, and that preceeds the
#  last period (the filename extension).
#  For instance, the tags of \code{path/to/foo,a.2,b.ext} are 
#  \code{a.2} and \code{b}.
# }
#
# @author
#
# \seealso{
#   @seemethod "getName".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "GenericDataFile", function(this, pattern=NULL, collapse=NULL, ...) {
  fullname <- getFullName(this, ...);

  # Data-set name is anything before the first comma
  name <- gsub("[,].*$", "", fullname);

  # Keep anything after the data-set name (and the separator).
  tags <- substring(fullname, nchar(name)+2);
  tags <- unlist(strsplit(tags, split=","));

  # Keep only those matching a regular expression?
  if (!is.null(pattern))
    tags <- grep(pattern, tags, value=TRUE);

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }
 
  if (length(tags) == 0)
    tags <- NULL;

  tags;
})

setMethodS3("hasTags", "GenericDataFile", function(this, tags, ...) {
  all(tags %in% getTags(this));
})

setMethodS3("hasTag", "GenericDataFile", function(this, tag, ...) {
  hasTags(this, tags=tag, ...);
})



###########################################################################/**
# @RdocMethod getFileType
#
# @title "Gets the file type of a file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character in lower case letters.
# }
#
# \details{
#   By default, this methods returns the filename extension, but subclasses
#   may override this.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFileType", "GenericDataFile", function(this, ...) {
  pattern <- "(.*)[.]([abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0-9]+)$";
  ext <- gsub(pattern, "\\2", this$.pathname);
  tolower(ext);
})


setMethodS3("isFile", "GenericDataFile", function(this, ...) {
  isFile(getPathname(this));
})


setMethodS3("getFileSize", "GenericDataFile", function(this, what=c("numeric", "units"), sep="", ...) {
  # Argument 'what':
  what <- match.arg(what);

  fileSize <-   file.info(this$.pathname)$size;
  if (what == "numeric")
    return(fileSize);

  if (is.na(fileSize))
    return(fileSize);
    
  units <- c("bytes", "kB", "MB", "GB", "TB");
  scale <- 1;
  for (kk in seq(along=units)) {
    unit <- units[kk];
    if (fileSize < 1000)
      break;
    fileSize <- fileSize/1024;
  }
  fileSize <- sprintf("%.2f%s%s", fileSize, sep, unit);

  fileSize;
})


setMethodS3("fromFile", "GenericDataFile", function(static, filename, path=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (.checkArgs) {
    # Argument 'filename' and 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
  } else {
    pathname <- filename;
  }


  # Get all known subclasses (bottom up)
  clazz <- Class$forName(class(static)[1]);
  knownSubclasses <- rev(getKnownSubclasses(clazz));
  for (className in knownSubclasses) {
    clazz <- Class$forName(className);

    # Try reading the file using the static fromFile() method of each class
    static <- getStaticInstance(clazz);
    tryCatch({
      res <- fromFile(static, filename=pathname, .checkArgs=FALSE);
      return(res);
    }, error = function(ex) {})
  }

  # If not "read" above, just create an instance as is.
  res <- newInstance(static, filename=pathname, ...);
  return(res);

##  throw("Could not read file '", pathname, "'.", 
##             "Tried using all known subclasses of ", class(static)[1], ": ", 
##                                     paste(knownSubclasses, collapse=", "));
}, static=TRUE)



setMethodS3("copyTo", "GenericDataFile", function(this, filename=getFilename(this), path=NULL, overwrite=FALSE, ..., verbose=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Assert that we're not trying to copy to itself
  if (identical(pathname, getPathname(this)))
    throw("Cannot copy file. Source and destination are identical: ", pathname);

  # Assert that file is not overwritten by mistake.
  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=!overwrite);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fail-safe copying
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Copying file");
  copyFile(getPathname(this), pathname, overwrite=overwrite, verbose=less(verbose, 10));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create object of the same class.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- newInstance(this, pathname);

  res;
}, protected=TRUE)



setMethodS3("renameTo", "GenericDataFile", function(this, filename=getFilename(this), path=NULL, ..., verbose=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Nothing to do?
  if (identical(pathname, getPathname(this)))
    return(this);

  # Assert that file is not overwritten by mistake.
  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=TRUE);

  srcPathname <- getPathname(this);

  verbose && enter(verbose, "Renaming ", class(this)[1], " pathname");
  verbose && cat(verbose, "Source: ", srcPathname);
  verbose && cat(verbose, "Destination: ", pathname);

  verbose && enter(verbose, "Renaming file");
  res <- file.rename(srcPathname, pathname);
  if (!res) {
    throw("Failed to rename file: ", srcPathname, " -> ", pathname);
  }
  verbose && exit(verbose);

  # Update GenericDataFile object
  this$.pathname <- pathname;

  verbose && exit(verbose);

  invisible(this);
}, protected=TRUE)



setMethodS3("getChecksum", "GenericDataFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Calculating checksum");
  pathname <- getPathname(this);
  checksum <- digest2(pathname, file=TRUE);
  verbose && exit(verbose);

  checksum;
})


setMethodS3("writeChecksum", "GenericDataFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  pathname <- getPathname(this);
  outPathname <- sprintf("%s.md5", pathname);

  verbose && enter(verbose, "Writing checksum");
  verbose && cat(verbose, "Pathname: ", outPathname);
  checksum <- getChecksum(this, verbose=less(verbose));
  cat(checksum, file=outPathname);
  verbose && exit(verbose);

  invisible(outPathname);
})



setMethodS3("readChecksum", "GenericDataFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  pathname <- getPathname(this);
  outPathname <- sprintf("%s.md5", pathname);

  verbose && enter(verbose, "Reading checksum");
  checksum <- readLines(outPathname, warn=FALSE);
  verbose && exit(verbose);

  checksum;
})


setMethodS3("compareChecksum", "GenericDataFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  pathname <- getPathname(this);
  outPathname <- sprintf("%s.md5", pathname);

  verbose && enter(verbose, "Comparing checksum");
  verbose && cat(verbose, "Pathname: ", outPathname);

  checksum <- getChecksum(this, verbose=less(verbose));
  if (isFile(outPathname)) {
    checksum2 <- readLines(outPathname, warn=FALSE);
  } else {
    checksum2 <- NA;
  }
  res <- identical(checksum, checksum2);

  verbose && cat(verbose, res);
  verbose && exit(verbose);

  res;
})


setMethodS3("validateChecksum", "GenericDataFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Validating checksum");
  pathname <- getPathname(this);
  res <- compareChecksum(this, ..., verbose=less(verbose));
  if (!res) {
    throw("The calculated checksum and the checksum store on file do not match: ", pathname);
  }
  verbose && exit(verbose);

  invisible(res);
})


setMethodS3("renameToUpperCaseExt", "GenericDataFile", function(static, pathname, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # isFileCase() is a case-sensitive isFile() for Windows
  isFileCase <- function(pathname, ...) {
    # Non-case sensitive check
    if (!isFile(pathname))
      return(FALSE);

    # There can still be a case-difference
    path <- dirname(pathname);
    filename <- basename(pathname);
    filenames <- list.files(path=path, all.files=TRUE);
    res <- grep(filename, filenames, fixed=TRUE);
    res <- (length(res) >= 1);
    res;
  }


  # Identify the filename extension
  ext <- gsub("^.*[.]", "", pathname);

  # No filename extension?  Do nothing.
  if (identical(ext, pathname))
    return(pathname);

  # Generate the pathname with lower- and upper-case extensions
  extL <- tolower(ext);
  pathnameL <- gsub(sprintf("[.]%s$", ext), sprintf(".%s", extL), pathname);

  extU <- toupper(ext);
  pathnameU <- gsub(sprintf("[.]%s$", ext), sprintf(".%s", extU), pathname);

  # Does a lower-case filename exist? If not, nothing to do.
  if (!isFileCase(pathnameL))
    return(pathnameU);

  # Can we rename it?
  if (identical(pathname, pathnameL) && isFileCase(pathnameU)) {
    throw("Cannot rename pathname to have upper-case filename extension, because such a file already exists: ", pathnameU);
  }

  # Try to rename the file
  res <- file.rename(pathnameL, pathnameU);
  if (res) {
    msg <- paste("Renamed file to have an upper-case filename extension:", pathname);
    warning(msg);
  } else {
    throw("Failed to rename file such that it gets an upper-case filename extension (try to rename the file manually): ", pathname);
  }

  pathnameU;
}, static=TRUE, protected=TRUE)


############################################################################
# HISTORY:
# 2008-05-11
# o Now static fromFile() always creates an instance.
# 2008-05-09
# o Moved private get/seLabel() from AffymetrixFile to GenericDataFile.
# o Moved the attributes features from AffymetrixFile to GenericDataFile.
# 2008-03-22
# o Added 'aliased' to getFullName().
# 2007-09-25
# o Added isFile() to test if the file exists or not.
# 2007-09-15
# o Added renameTo().
# o Now copyTo() utilizes fileCopy(), which is fail safe.
# 2007-09-14
# o Extracted GenericDataFile from AffymetrixFile.
# 2007-09-13
# o Added missing setAttributesByTags().
# 2007-08-09
# o Added static renameToUpperCaseExt().
# 2007-03-20
# o Added getAlias() and setAlias().  Note, getName() etc are still
#   unaffected by these.
# 2007-03-05
# o Added setAttributesByTags(), which now also tries to coerce values.
# o Added support for (in-memory) attributes.
# 2007-02-07
# o Added getChecksum(), writeChecksum(), readChecksum(), and 
#   compareChecksum() and validateChecksum(). I did this because I noticed 
#   by chance that some of my CEL files transferred via an external HDD got
#   corrupt probe signals.
# 2007-01-14
# o Added a test for "unknown" (=unused) arguments to constructor.
# 2007-01-07
# o Added hasTags() and hasTag().
# 2006-11-02
# o Added getFullName(), getTags() and redefined getName().
# 2006-09-15
# o Added stextSize().
# 2006-08-27
# o Added stextLabel() and stextLabels(). stext is for "side text", cf. 
#   mtext for "margin text". stext() is slightly more convenient than mtext
#   when it comes to different font sizes.
# o Added copyTo().
# 2006-08-14
# o Added abstract fromFile().
# 2006-08-11
# o Created from AffymetrixDataFile in order to represent CDF files too.
############################################################################
