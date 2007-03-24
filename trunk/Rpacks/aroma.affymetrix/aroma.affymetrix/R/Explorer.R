###########################################################################/**
# @RdocClass Explorer
#
# @title "The Explorer class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{tags}{A @character @vector of tags to be added to the output path.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
# }
#*/###########################################################################
setConstructorS3("Explorer", function(tags="*", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
  }


  extend(Object(), "Explorer",
    .alias = NULL,
    .tags = tags,
    .arrays = NULL
  )
})



setMethodS3("as.character", "Explorer", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Main path: %s", getMainPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)



###########################################################################/**
# @RdocMethod getArrays
#
# @title "Gets the names of the arrays"
#
# \description{
#  @get "title" in the explorer.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getArrays", "Explorer", function(this, ...) {
  arrays <- this$.arrays;
  if (is.null(arrays)) {
    arrays <- getArraysOfInput(this);
  }
  arrays;
})



###########################################################################/**
# @RdocMethod setArrays
#
# @title "Sets the arrays"
#
# \description{
#  @get "title" to be processed by the explorer.
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{A @character (or @integer) @vector of arrays to be 
#      considered. If @NULL, all arrays of the data set are considered.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setArrays", "Explorer", abstract=TRUE);



###########################################################################/**
# @RdocMethod nbrOfArrays
#
# @title "Gets the total number of arrays"
#
# \description{
#  @get "title" considered by the explorer.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfArrays", "Explorer", function(this, ...) {
  length(getArrays(this));
})



###########################################################################/**
# @RdocMethod getAlias
#
# @title "Gets the alias of the output set"
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
setMethodS3("getAlias", "Explorer", function(this, ...) {
  this$.alias;
})



###########################################################################/**
# @RdocMethod setAlias
#
# @title "Sets the alias of the output set"
#
# \description{
#   @get "title".
#   If specified, the alias overrides the data set name, which is used by 
#   default.
# }
#
# @synopsis
#
# \arguments{
#  \item{alias}{A @character string for the new alias of the output set.
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
setMethodS3("setAlias", "Explorer", function(this, alias=NULL, ...) {
  # Argument 'alias':
  if (!is.null(alias)) {
    alias <- Arguments$getFilename(alias);  # Valid filename?

    # Assert that no commas are used.
    if (regexpr("[,]", alias) != -1) {
      throw("Output-set aliases (names) must not contain commas: ", alias);
    }
  }

  this$.alias <- alias;
})


###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the explorer"
#
# \description{
#  @get "title", which is the same as the name of the data set.
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
# \details{
#  If a name alias has not been set explicitly, the name of the data set will
#  used.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "Explorer", function(this, ...) {
  name <- getAlias(this);
  if (is.null(name)) {
    name <- getNameOfInput(this);
  }
  name;
})



###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the explorer"
#
# \description{
#  @get "title", which are the tags of the data set plus additional tags.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "Explorer", function(this, collapse=NULL, ...) {
  tags <- getTagsOfInput(this, ...);

  tags <- c(tags, this$.tags);
  tags <- unique(tags);

  # Update asterisk tags
  tags[tags == "*"] <- getAsteriskTags(this);

  # Keep non-empty tags
  tags <- tags[nchar(tags) > 0];

  tags <- unique(tags);

  tags <- paste(tags, collapse=collapse);
  if (length(tags) == 0)
    tags <- NULL;

  tags;
})

setMethodS3("getAsteriskTags", "Explorer", function(this, ...) {
  "";
})

setMethodS3("getTagsOfInput", "Explorer", function(this, ...) {
  "";
})

setMethodS3("getNameOfInput", "Explorer", abstract=TRUE, protected=TRUE);



setMethodS3("getFullName", "Explorer", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})



setMethodS3("getRootPath", "Explorer", function(this, ...) {
  "reports";
}, private=TRUE)


setMethodS3("getMainPath", "Explorer", function(this, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  name <- getName(this);

  # Tags
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  if (nchar(tags) == 0) {
    tags <- "raw";  # Default
  }

  # The full path
  path <- filePath(rootPath, name, tags, expandLinks="any");
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
}, protected=TRUE)

setMethodS3("getPath", "Explorer", abstract=TRUE);


setMethodS3("getTemplatePath", "Explorer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating template files for ChromosomeExplorer");
  # Search for template files
  rootPath <- getRootPath(this);
  path <- filePath(rootPath, "templates", expandLinks="any");
  if (!isDirectory(path)) {
    path <- system.file("reports", "templates", package="aroma.affymetrix");
  }
  verbose && exit(verbose);

  path;
}, protected=TRUE)


setMethodS3("getIncludePath", "Explorer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating include files for ChromosomeExplorer");
  # Search for include files
  path <- system.file("reports", "includes", package="aroma.affymetrix");
  verbose && exit(verbose);

  path;
}, protected=TRUE)


setMethodS3("addIncludes", "Explorer", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Setting up ", class(this)[1], " report files");

  destPath <- filePath(getRootPath(this), "includes");
  if (force || !isDirectory(destPath)) {
    verbose && enter(verbose, "Copying template files");
    srcPath <- getIncludePath(this);
    verbose && cat(verbose, "Source path: ", srcPath);
    verbose && cat(verbose, "Destination path: ", destPath);
    pathnames <- copyDirectory(from=srcPath, to=destPath, recursive=TRUE);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);
})


setMethodS3("addIndexFile", "Explorer", function(this, filename=sprintf("%s.html", class(this)[1]), ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  srcPath <- getTemplatePath(this);
  srcPathname <- filePath(srcPath, "html", class(this)[1], filename);
  outPathname <- filePath(getMainPath(this), filename);

  if (force || !isFile(outPathname)) {
    verbose && enter(verbose, "Copying ", filename);
    verbose && cat(verbose, "Source pathname: ", srcPathname);
    verbose && cat(verbose, "Destination pathname: ", outPathname);
    if (!isFile(srcPathname))
      throw("File not found: ", srcPathname);
    file.copy(srcPathname, outPathname, overwrite=TRUE);
    verbose && exit(verbose);
  }
}, protected=TRUE)



setMethodS3("setup", "Explorer", abstract=TRUE);



###########################################################################/**
# @RdocMethod process
#
# @title "Generates image files, scripts and dynamic pages for the explorer"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{A @vector of arrays specifying which arrays to
#    be considered.  If @NULL, all are processed.}
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "Explorer", abstract=TRUE);


###########################################################################/**
# @RdocMethod display
#
# @title "Displays the explorer in the default browser"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("display", "Explorer", function(this, filename=sprintf("%s.html", class(this)[1]), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Opening ", class(this)[1]);

  # The path to the explorer HTML document
  path <- getMainPath(this);
  pathname <- filePath(path, filename, expandLinks="any");

  # Just in case, is setup needed?
  if (!isFile(pathname)) {
    setup(this, verbose=less(verbose));
    if (!isFile(pathname)) {
      throw("Cannot open ", class(this)[1], ". No such file: ", pathname);
    }
  }

  pathname <- getAbsolutePath(pathname);
  pathname <- chartr("/", "\\", pathname);

  verbose && cat(verbose, "Pathname: ", pathname);
  res <- browseURL(pathname, ...);

  verbose && exit(verbose);

  invisible(res);
})


##############################################################################
# HISTORY:
# 2007-03-24
# o BUG FIX: getPath() created the root path before trying to expand
#   Windows shortcuts.
# 2007-03-19
# o Created from ChromosomeExplorer.R.
##############################################################################
