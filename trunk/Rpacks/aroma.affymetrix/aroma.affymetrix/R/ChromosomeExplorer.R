###########################################################################/**
# @RdocClass ChromosomeExplorer
#
# @title "The ChromosomeExplorer class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{model}{A @see "GladModel" object.}
#   \item{tags}{A @character @vector of tags to be added to the output path.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Generating PNG images}{
#   In order to get better looking graphs, but also to be able to generate
#   bitmap images on systems without direct bitmap support, which is the case
#   when running R in batch mode or on Unix without X11 support, images are
#   created using the @see "R.utils::png2" device (a wrapper for 
#   \code{bitmap()} immitating \code{png()}).  The \code{png()} is only
#   used if \code{png2()}, which requires Ghostscript, does not.
#   Note, when images are created using \code{png2()}, the images does
#   not appear immediately, although the function call is completed,
#   so be patient.
# }
#
# @author
# 
# \seealso{
#  @see "GladModel".
# }
#*/###########################################################################
setConstructorS3("ChromosomeExplorer", function(model=NULL, tags="*", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  if (!is.null(model)) {
    if (!inherits(model, "GladModel")) {
      throw("Argument 'model' is not a 'GladModel': ", class(model));
    }
  }

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    tags[tags == "*"] <- "";
  }


  extend(Object(), "ChromosomeExplorer",
    .model = model
  )
})


setMethodS3("as.character", "ChromosomeExplorer", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getModel", "ChromosomeExplorer", function(this, ...) {
  this$.model;
})

setMethodS3("nbrOfArrays", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  nbrOfArrays(model);
})

setMethodS3("getArrays", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  getArrays(model, ...);
})


setMethodS3("getName", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  getName(model, ...);
})



setMethodS3("getTags", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  tags <- getTags(model);
  tags <- c(tags, this$.tags);
  tags <- unique(tags);
  tags;
})


setMethodS3("getFullName", "ChromosomeExplorer", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "ChromosomeExplorer", function(this, ...) {
  "ce";
}, private=TRUE)


setMethodS3("getPath", "ChromosomeExplorer", function(this, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  name <- getName(this);

  # Tags
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");

  # Chip type    
  model <- getModel(this);
  chipType <- getChipType(model);

  # The full path
  path <- filePath(rootPath, name, tags, chipType, expandLinks="any");
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
})


###########################################################################/**
# @RdocMethod getChromosomes
#
# @title "Gets the chromosomes available"
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
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChromosomes", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  getChromosomes(model);
})


setMethodS3("writeGraphs", "ChromosomeExplorer", function(x, ...) {
  # To please R CMD check.
  this <- x;

  path <- getPath(this);
  model <- getModel(this);
  plot(model, path=path, imageFormat="png", ...);

  invisible(path);
}, private=TRUE)

setMethodS3("writeRegions", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
#  writeRegions(model, ...);
}, private=TRUE)


setMethodS3("process", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  fit(model, ...);
  writeGraphs(this, ...);
  writeRegions(this, ...);
})


##############################################################################
# HISTORY:
# 2007-01-04
# o Created.
##############################################################################
