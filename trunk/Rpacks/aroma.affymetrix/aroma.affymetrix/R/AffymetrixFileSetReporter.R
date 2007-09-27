###########################################################################/**
# @RdocClass AffymetrixFileSetReporter
#
# @title "The AffymetrixFileSetReporter class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{set}{An @see "AffymetrixFileSet" object.}
#   \item{tags}{A @character @vector of tags to be added to the output path.}
#   \item{...}{Not used.}
#   \item{.setClass}{The name of the class of the input set.}
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
setConstructorS3("AffymetrixFileSetReporter", function(set=NULL, tags="*", ..., .setClass="AffymetrixFileSet") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'set':
  if (is.null(set)) {
  } else if (!inherits(set, .setClass)) {
    throw("Argument 'set' is not a ", .setClass, ": ", class(set)[1]);
  }

  extend(Object(...), "AffymetrixFileSetReporter",
    .alias = NULL,
    .set = set
  )
})


setMethodS3("as.character", "AffymetrixFileSetReporter", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getFileSet", "AffymetrixFileSetReporter", function(this, ...) {
  this$.set;
}, protected=TRUE)



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
setMethodS3("getAlias", "AffymetrixFileSetReporter", function(this, ...) {
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
setMethodS3("setAlias", "AffymetrixFileSetReporter", function(this, alias=NULL, ...) {
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
setMethodS3("getName", "AffymetrixFileSetReporter", function(this, ...) {
  name <- getAlias(this);
  if (is.null(name)) {
    set <- getFileSet(this);
    name <- getName(set);
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
setMethodS3("getTags", "AffymetrixFileSetReporter", function(this, collapse=NULL, ...) {
  set <- getFileSet(this);
  tags <- getTags(set);

  tags <- c(tags, this$.tags);
  tags <- locallyUnique(tags);

  # Update asterisk tags
  tags[tags == "*"] <- getAsteriskTags(this);

  # Keep non-empty tags
  tags <- tags[nchar(tags) > 0];

  tags <- locallyUnique(tags);

  tags <- paste(tags, collapse=collapse);
  if (length(tags) == 0)
    tags <- NULL;

  tags;
})

setMethodS3("getAsteriskTags", "AffymetrixFileSetReporter", function(this, ...) {
  "";
})


setMethodS3("getFullName", "AffymetrixFileSetReporter", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getReportSet", "AffymetrixFileSetReporter", abstract=TRUE, protected=TRUE);


setMethodS3("getRootPath", "AffymetrixFileSetReporter", function(this, ...) {
  "reports";
}, private=TRUE)


setMethodS3("getMainPath", "AffymetrixFileSetReporter", function(this, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  name <- getName(this);

  # Tags
  tags <- getTags(this, collapse=",");
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


setMethodS3("getPath", "AffymetrixFileSetReporter", abstract=TRUE);


setMethodS3("setup", "AffymetrixFileSetReporter", abstract=TRUE);


###########################################################################/**
# @RdocMethod process
#
# @title "Generates report"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{See subclasses.}
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
setMethodS3("process", "AffymetrixFileSetReporter", abstract=TRUE);



##############################################################################
# HISTORY:
# 2007-03-24
# o BUG FIX: getPath() created the root path before trying to expand
#   Windows shortcuts.
# 2007-03-19
# o Created from Explorer.R.
##############################################################################
