###########################################################################/**
# @RdocClass AffymetrixFileSet
#
# @title "The AffymetrixFileSet class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixFileSet object represents a set of @see "AffymetrixFile"s
#  with \emph{identical} chip types.
# }
# 
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "AffymetrixFile":s.}
#   \item{tags}{A @character @vector of tags to be used for this file set.
#      The string \code{"*"} indicates that it should be replaced by the
#      tags part of the file set pathname.}
#   \item{alias}{A @character string specifying a name alias overriding the
#      name inferred from the pathname.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("AffymetrixFileSet", function(files=NULL, tags="*", alias=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'files':
  if (is.null(files)) {
  } else if (is.list(files)) {
    lapply(files, FUN=function(df) {
      if (!inherits(df, "AffymetrixFile"))
        throw("Argument 'files' contains a non-AffymetrixFile object: ", class(df));
    })
  } else if (inherits(files, "AffymetrixFileSet")) {
    return(as.AffymetrixFileSet(files));
  } else {
    throw("Argument 'files' is of unknown type: ", mode(files));
  }

  # Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  }


  this <- extend(Object(), "AffymetrixFileSet",
    files = as.list(files),
    .alias = NULL,
    .tags = NULL
  );

  setTags(this, tags);
  setAlias(this, alias);

  this;
})


###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the Affymetrix file set"
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
setMethodS3("as.character", "AffymetrixFileSet", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Name: %s", getName(this)));
  tags <- getTags(this, collapse=",");
  if (!is.null(tags)) {
    s <- paste(s, " Tags: ", tags, ".", sep="");
  }
  s <- c(s, sprintf("Full name: %s", getFullName(this)));
  n <- nbrOfFiles(this);
  s <- c(s, sprintf("Number of files: %d", n));
  names <- getNames(this);
  if (n >= 5)
    names <- c(names[1:2], "...", names[n]);
  names <- paste(names, collapse=", ");
  s <- c(s, sprintf("Names: %s", names));
  s <- c(s, sprintf("Path (to the first file): %s", getPath(this)));
  s <- c(s, sprintf("Total file size: %.2fMB", getFileSize(this)/1024^2));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)




setMethodS3("clone", "AffymetrixFileSet", function(this, clear=TRUE, ...) {
  # Clone itself
  object <- NextMethod("clone", this, ...);

  # Clone each file object
  files <- as.list(object);
  for (kk in seq(along=files)) {
    files[[kk]] <- clone(files[[kk]], clear=TRUE);
  }
  object$files <- files;

  # Clear the cached fields?
  if (clear)
    clearCache(object);

  object;
}, private=TRUE)



###########################################################################/**
# @RdocMethod getFullName
#
# @title "Gets the full name of the file set"
#
# \description{
#   @get "title", that is the name of the directory without parent directories.
# }
#
# @synopsis
#
# \arguments{
#  \item{parent}{The number of generations up in the directory tree the
#    directory name should be retrieved.  By default the current directory
#    is used.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character.
# }
#
# \details{
#  By default, the full name of a file set is the name of the directory 
#  containing all the files, e.g. the name of file set \code{path/to,a,b/*} 
#  is \code{to,a,b}.
#  Argument \code{parent=1} specifies that the parent directory should be
#  used, and so on.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFullName", "AffymetrixFileSet", function(this, parent=1, ...) {
  parent <- Arguments$getInteger(parent, range=c(0,32));

  # The name of a file set is inferred from the pathname of the directory
  # of the set assuming path/to/<name>/<something>/<"chip type">/

  # Get the path of this file set
  path <- getPath(this);

  while (parent > 0) {
    # path/to/<name>/<something>
    path <- dirname(path);
    parent <- parent - 1;
  }

  # <name>
  name <- basename(path);

  name;
})


###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the file set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional arguments passed to @seemethod "getFullName", if
#    the name is inferred from the pathname.}
# }
#
# \value{
#   Returns a @character.
# }
#
# \details{
#  If a name alias has not been set explicitly, the name is inferred from 
#  the pathname of the file set.
#  The \emph{name} of a file set is the part of the directory name that 
#  preceeds the first comma, if any.
#  For instance, the name of the file set named \code{foo,a,b} is \code{foo}.
# }
#
# @author
#
# \seealso{
#   @seemethod "getAlias".
#   @seemethod "getFullName".
#   @seemethod "getTags".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "AffymetrixFileSet", function(this, ...) {
  # Use name alias, if specified
  name <- getAlias(this);
  if (is.null(name)) {
    # ...otherwise, infer from pathname.
    name <- getFullName(this, ...);

    # Keep anything before the first comma
    name <- gsub("[,].*$", "", name);
  }
  
  name;
})




###########################################################################/**
# @RdocMethod getAlias
#
# @title "Gets the alias of the file set"
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
#   Returns a @character, or @NULL if no alias is set.
# }
#
# @author
#
# \seealso{
#   @seemethod "setAlias".
#   @seemethod "getName".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAlias", "AffymetrixFileSet", function(this, ...) {
  this$.alias;
})



###########################################################################/**
# @RdocMethod setAlias
#
# @title "Sets the alias of the file set"
#
# \description{
#   @get "title".
#   If specified, the alias overrides the name inferred from the pathname
#   of the file set.  This can be used in order to use another name of the
#   output data set than the input data set of many transforms and models.
# }
#
# @synopsis
#
# \arguments{
#  \item{alias}{A @character string for the new alias of the file set.
#   The alias must consists of valid filename characters, and must not
#   contain commas, which are used to separate tags.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns nothing.
# }
#
# \details{
# }
#
# @author
#
# \seealso{
#   @seemethod "getAlias".
#   @seemethod "getName".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setAlias", "AffymetrixFileSet", function(this, alias=NULL, ...) {
  # Argument 'alias':
  if (!is.null(alias)) {
    alias <- Arguments$getFilename(alias);  # Valid filename?

    # Assert that no commas are used.
    if (regexpr("[,]", alias) != -1) {
      throw("File-set aliases (names) must not contain commas: ", alias);
    }
  }

  this$.alias <- alias;
})




###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the file set"
#
# \description{
#   @get "title".
#   Any tag equals \code{"*"} is replaced by the comma separated tags part of
#   the file-set pathname.
# }
#
# @synopsis
#
# \arguments{
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
#  The \emph{tags} of a file set are the comma separated parts of the
#  filename that follows the the first comma, if any, and that preceeds the
#  last period (the filename extension).
#  For instance, the tags of \code{path/to/foo,a.2,b.ext} are 
#  \code{a.2} and \code{b}.
# }
#
# @author
#
# \seealso{
#   @seemethod "setTags".
#   @seemethod "getName".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "AffymetrixFileSet", function(this, collapse=NULL, ...) {
  tags <- this$.tags;

  if ("*" %in% tags) {
    name <- getFullName(this, ...);

    # File-set name is anything before the first comma
    fsName <- gsub("[,].*$", "", name);

    # Keep anything after the file-set name (and the separator).
    name <- substring(name, nchar(fsName)+2);
  
    filenameTags <- strsplit(name, split=",")[[1]];

    pos <- which("*" == tags);
    tags <- tags[-pos];
    if (length(filenameTags) > 0) {
      if (length(tags) == 0) {
        tags <- filenameTags;
      } else {
        tags <- R.utils::insert.default(tags, pos[1], filenameTags); 
      }
    }
  }

  tags <- paste(tags, collapse=collapse);
  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("hasTags", "AffymetrixFileSet", function(this, tags, ...) {
  all(tags %in% getTags(this));
})


setMethodS3("hasTag", "AffymetrixFileSet", function(this, tag, ...) {
  hasTags(this, tags=tag, ...);
})



###########################################################################/**
# @RdocMethod setTags
#
# @title "Sets the tags of the file set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{tags}{A @character @vector of tags.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns nothing.
# }
#
# \details{
#   See @seemethod "getTags" for so called \emph{special tags}.
# }
#
# @author
#
# \seealso{
#   @seemethod "getTags".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setTags", "AffymetrixFileSet", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
  }
  
  this$.tags <- tags;
})


setMethodS3("getFileSize", "AffymetrixFileSet", function(this, ...) {
  if (is.null(fileSize <- this$.fileSize)) {
    fileSize <- sum(unlist(lapply(this, FUN=getFileSize), use.names=FALSE));
    this$.fileSize <-  fileSize;
  }
  fileSize;
})



###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path (directory) of the file set"
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
setMethodS3("getPath", "AffymetrixFileSet", function(this, ...) {
  getPath(this$files[[1]]);
})



###########################################################################/**
# @RdocMethod length
# @aliasmethod nbrOfFiles
#
# @title "Gets the number of files in the file set"
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
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("length", "AffymetrixFileSet", function(this, ...) {
  length(this$files);
}, private=TRUE)

setMethodS3("nbrOfFiles", "AffymetrixFileSet", function(this, ...) {
  length(this, ...);
})



###########################################################################/**
# @RdocMethod lapply
#
# @title "Applies a function to each of the Affymetrix file files"
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
#   Returns a named @list where the names are the names of the date files
#   according to @seemethod "getNames".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("lapply", "AffymetrixFileSet", function(this, ...) {
  res <- lapply(this$files, ...);
  names(res) <- unlist(lapply(as.list(this), FUN=getFullName));
  res;
})

setMethodS3("sapply", "AffymetrixFileSet", function(this, ...) {
  res <- sapply(this$files, ...);
  names(res) <- unlist(lapply(as.list(this), FUN=getFullName));
  res;
})


setMethodS3("reorder", "AffymetrixFileSet", function(x, order, ...) {
  # To please R CMD check
  this <- x;

  # Argument 'order':
  if (is.character(order)) {
    # Assume 'order' contains names
    names <- getNames(this);

    # Identify special tags and remove them
    idx <- (order == "*");
    if (sum(idx) > 1)
      throw("Argument 'order' contains more than one asterix.");
    pos <- match(order[!idx], names);
    if (any(is.na(pos))) {
      bad <- order[!idx][is.na(pos)];
      throw("Argument 'order' contains unknown sample names: ", 
                                                 paste(bad, collapse=", "));
    }
    if (sum(idx) == 0) {
      order <- pos;
    } else {
      order <- as.list(order);
      order[!idx] <- names[pos];
      order[[which(idx)]] <- setdiff(names, names[pos]);
      order <- unlist(order, use.names=FALSE);
      order <- match(order, names);
    }
  }

  order <- Arguments$getIndices(order, range=c(1, nbrOfFiles(this)));
  if (any(duplicated(order))) {
    bad <- order[duplicated(order)];
    throw("Argument 'order' contains duplicates: ", 
                                                 paste(bad, collapse=", "));
  }

  this$files <- this$files[order];
  invisible(this);
}, private=TRUE)


###########################################################################/**
# @RdocMethod getNames
#
# @title "Gets the names of the files in the file set"
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
#   Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getNames", "AffymetrixFileSet", function(this, ...) {
  res <- unlist(lapply(this, FUN=getName));
  unname(res);
})

setMethodS3("getFullNames", "AffymetrixFileSet", function(this, ...) {
  res <- unlist(lapply(this, FUN=getFullName));
  unname(res);
})



###########################################################################/**
# @RdocMethod getPathnames
#
# @title "Gets the pathnames of the files in the file set"
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
#   Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPathnames", "AffymetrixFileSet", function(this, ...) {
  unlist(lapply(this, FUN=getPathname))
})



###########################################################################/**
# @RdocMethod seq
#
# @title "Gets an vector of file indices"
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
#   Returns an @integer @vector in [1,N] where N is the number of files,
#   or an empty vector if the file set is empty.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("seq", "AffymetrixFileSet", function(this, ...) {
  seq(length=nbrOfFiles(this));
})


###########################################################################/**
# @RdocMethod as.list
#
# @title "Returns the files of the Affymetrix file set"
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
#  Returns a @list of @see "AffymetrixFile"s.
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
setMethodS3("as.list", "AffymetrixFileSet", function(x, ...) {
  # To please R CMD check.
  this <- x;

  this$files;
})

setMethodS3("getFile", "AffymetrixFileSet", function(this, idx, ...) {
  this$files[[idx]];
})

setMethodS3("getFiles", "AffymetrixFileSet", function(this, idxs=NULL, ...) {
  if (is.null(idxs)) {
    this$files;
  } else {
    this$files[idxs];
  }
}, private=TRUE)


setMethodS3("appendFiles", "AffymetrixFileSet", function(this, files, clone=TRUE, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Appending ", length(files), " files");
  # Validate classes
  verbose && enter(verbose, "Validating file classes");
  className <- class(this$files[[1]])[1];
  isValid <- unlist(lapply(files, FUN=inherits, className));
  if (!all(isValid)) {
    throw("Some of the elements in argument 'files' are not '", 
      className, "'");
  }
  verbose && exit(verbose);

  # Clone file objects?
  if (clone) {
    verbose && enter(verbose, "Cloning files");
    files <- lapply(files, FUN=function(file) clone(file));    
    verbose && exit(verbose);
  }

  # Append
  this$files <- append(this$files, files);

  verbose && exit(verbose);

  invisible(this);
})


setMethodS3("append", "AffymetrixFileSet", function(x, values, ...) {
  # To please R CMD check
  this <- x;
  other <- values;

  if (!inherits(other, class(this)[1])) {
    throw("Argument 'other' is not an ", class(this)[1], " object: ", 
                                                      class(other)[1]);
  }

  appendFiles(this, getFiles(other), ...);
})



###########################################################################/**
# @RdocMethod extract
#
# @title "Extract a subset of the file set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{files}{An @integer or a @logical @vector indicating which data files
#    to be extracted.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @see "AffymetrixFileSet" (or a subclass) object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extract", "AffymetrixFileSet", function(this, files, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.logical(files))
    files <- which(files);
  files <- Arguments$getIndices(files, range=range(seq(this)));

  res <- clone(this);
  res$files <- this$files[files];
  clearCache(res);  # Some cached values are incorrect now.

  res;
})


###########################################################################/**
# @RdocMethod as.AffymetrixFileSet
# @alias as.AffymetrixFileSet.list
# @alias as.AffymetrixFileSet.default
#
# @title "Coerce an object to an AffymetrixFileSet object"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Other arguments passed to @see "base::list.files".}
# }
#
# \value{
#   Returns an @see "AffymetrixFileSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("as.AffymetrixFileSet", "AffymetrixFileSet", function(object, ...) {
  object;
})

setMethodS3("as.AffymetrixFileSet", "list", function(object, ...) {
  AffymetrixFileSet(object, ...);
})

setMethodS3("as.AffymetrixFileSet", "default", function(object, ...) {
  throw("Cannot coerce object to an AffymetrixFileSet object: ", mode(object));
})


setMethodS3("clearCache", "AffymetrixFileSet", function(this, ...) {
  # Clear the cache of all files
  lapply(this, clearCache);

  # Clear the cache of the CDF object
  clearCache(getCdf(this));

  # Then for this object
  NextMethod("clearCache", this);
})




###########################################################################/**
# @RdocMethod fromFiles
#
# @title "Defines an AffymetrixFileSet object by searching for Affymetrix files"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{path}{The path where to search for Affymetrix files.}
#  \item{pattern}{The filename pattern for match files. 
#     If @NULL, filename extensions corresponding to known subclasses
#     of the abstract @see "AffymetrixFile" class are search for.}
#  \item{recursive}{If @TRUE, subdirectories are search recursively,
#     otherwise not.}
#  \item{...}{Optional arguments passed to the constructor of the
#     static (calling) class.}
# }
#
# \value{
#   Returns an @see "AffymetrixFileSet" object.
# }
#
# \section{Reserved filenames}{
#   Note that files with names starting with a period \code{.} are not 
#   searched for.  The reason for this is that such files are reserved for
#   internal use of this package.  For instance, the package store average
#   signals across CEL files in a file named as \code{.average<something>.CEL}
#   in the same directory as the CEL files, and when such a directory is
#   scanned we do not want such files to be interpreted as data.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromFiles", "AffymetrixFileSet", function(static, path=NULL, pattern=NULL, recursive=FALSE, fileClass="AffymetrixFile", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Argument 'pattern':
  if (!is.null(pattern))
    pattern <- Arguments$getRegularExpression(pattern);

  # Argument 'fileClass':
  clazz <- Class$forName(fileClass);
  dfStatic <- getStaticInstance(clazz);
  if (!inherits(dfStatic, "AffymetrixFile"))
    throw("Argument 'fileClass' is not refering to a AffymetrixFile class: ", fileClass);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Defining an ", class(static)[1], " object from files");
  verbose && cat(verbose, "Path: ", path);
  verbose && cat(verbose, "Pattern: ", pattern);
  verbose && cat(verbose, "File class: ", class(dfStatic)[1]);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create AffymetrixFile object from the matching files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Scan for Affymetrix files
  verbose && enter(verbose, "Scanning directory for files");
  pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE, 
                                   all.files=FALSE, recursive=recursive);
  verbose && printf(verbose, "Found %d files.\n", length(pathnames));
  if (length(pathnames) == 0)
    throw("No files found: ", path);
  verbose && exit(verbose);

  # Sort files in lexicographic order
  pathnames <- sort(pathnames);

  verbose && enter(verbose, "Defining ", length(pathnames), " files");
  files <- list();
  for (kk in seq(along=pathnames)) {
    if (as.logical(verbose)) cat(kk, ", ", sep="");
    df <- fromFile(dfStatic, pathnames[kk], .checkArgs=FALSE, verbose=less(verbose));
    files[[kk]] <- df;
    if (kk == 1) {
      clazz <- Class$forName(class(df)[1]);
      dfStatic <- getStaticInstance(clazz);
    }
  }
  if (as.logical(verbose)) cat("\n");
  verbose && exit(verbose);

  # Create the file set object
  set <- newInstance(static, files, ...);

  verbose && exit(verbose);

  set;
}, static=TRUE)



############################################################################
# HISTORY:
# 2007-02-15
# o Added getFullNames().
# 2007-02-07
# o Added sapply().
# 2007-01-14
# o Added a test for "unknown" (=unused) arguments to constructor.
# 2006-01-07
# o Added hasTags() and hasTag().
# o Now arguments '...' in fromFiles() are passed to the constructor of
#   the static class, i.e. via newInstance(static, ...).  Requested by KS.
# o Added argument 'alias' to constructor.
# o Added getAlias() and setAlias(), where the latter replaces setName().
# 2006-12-02
# o Added reorder().
# 2006-12-01
# o Now lapply() add data file names to returned list.
# 2006-11-20
# o Added support to override name of file set.
# o Added support for optional tags.
# 2006-11-02
# o Added getFullName(), getTags() and redefined getName().
# 2006-10-22
# o Now 'recursive' of fromFiles() defaults to FALSE.
# o Added getFiles() again.
# 2006-09-11
# o Added getPathnames().
# 2006-08-27
# o Added getFile() and getFiles().
# o Made filenames starting with a period reserved for internal use.
# 2006-08-26
# o Now getName() of a file set is inferred from the pathname:
#     path/to/<name>/chip_files/<"chip type">/
# 2006-08-21
# o Renamed 'array' to 'file'.
# o Extracted from AffymetrixCelSet.R.
# 2006-08-11
# o Added clearCache() which also clears the cache of all data file object.
# 2006-05-16
# o Redefined "[" to extract arrays.
# 2006-04-13
# o Added Rdoc comments for all methods.
# 2006-04-09
# o Now the read map is loaded automatically when fromFiles() used.
# 2006-03-30
# o Updated to new aroma.apd.
# 2006-03-18
# o Added argument 'subset' to calcAvgCellSignals() & normalizeQuantile().
# 2006-03-15
# o Now nbrOfCells() returns the number of cells for the first file only.
# o Now the fromFiles(static, ...) creates an object of the same class as 
#   the static object.
# 2006-03-04
# o Added mapping functions.
# o Added writeApd().
# 2006-03-03
# o Added lapply().
# 2006-03-02
# o Updated to deal with AffymetrixDataFile object instead of CEL directly.
# 2006-02-21
# o Letting readCelUnits() transform signals improves speed substantially.
# o Making use of new multi-array readCelUnits().
# 2006-02-20
# o Created.
############################################################################
