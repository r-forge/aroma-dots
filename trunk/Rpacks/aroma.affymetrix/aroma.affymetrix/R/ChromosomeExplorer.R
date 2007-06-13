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
setConstructorS3("ChromosomeExplorer", function(model=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  if (!is.null(model)) {
    if (!inherits(model, "GladModel")) {
      throw("Argument 'model' is not a 'GladModel': ", class(model)[1]);
    }
  }

  extend(Explorer(...), "ChromosomeExplorer",
    .model = model,
    .arrays = NULL,
    .plotCytoband = TRUE
  )
})


setMethodS3("as.character", "ChromosomeExplorer", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, paste("Number of arrays:", nbrOfArrays(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("setCytoband", "ChromosomeExplorer", function(this, status=TRUE, ...) {
  # Argument 'status':
  status <- Arguments$getLogical(status);

  this$.plotCytoband <- status;
})


###########################################################################/**
# @RdocMethod getModel
#
# @title "Gets the model"
#
# \description{
#  @get "title" for which the explorer is displaying it results.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "GladModel".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getModel", "ChromosomeExplorer", function(this, ...) {
  this$.model;
})


setMethodS3("getNames", "ChromosomeExplorer", function(this, ...) {
  getArrays(this, ...);
})

setMethodS3("getArraysOfInput", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  getNames(model, ...);
}, protected=TRUE)


setMethodS3("getFullNames", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  arrays <- getArrays(this);
  idx <- match(arrays, getNames(model));
  fullnames <- getFullNames(model, arrays=idx);
  fullnames;
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
#   \item{arrays}{A @character (or @integer) @vector of arrays used in the
#      model.  If @NULL, all arrays of the model are considered.}
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
setMethodS3("setArrays", "ChromosomeExplorer", function(this, arrays=NULL, ...) {
  # Argument 'arrays':
  if (!is.null(arrays)) {
    setTuple <- getSetTuple(this);
    arrays <- indexOfArrays(setTuple, arrays);
    arrays <- getArrays(setTuple)[arrays];
  }

  this$.arrays <- arrays;
})

setMethodS3("getSetTuple", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  getSetTuple(model);
}, protected=TRUE)


setMethodS3("getNameOfInput", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  getName(model, ...);
}, protected=TRUE)


setMethodS3("getTagsOfInput", "ChromosomeExplorer", function(this, ...) {
  model <- getModel(this);
  getTags(model);
}, protected=TRUE)


setMethodS3("getPath", "ChromosomeExplorer", function(this, ...) {
  mainPath <- getMainPath(this);

  # Chip type    
  model <- getModel(this);
  chipType <- getChipType(model);

  # Image set
  set <- "glad";

  # The full path
  path <- filePath(mainPath, chipType, set, expandLinks="any");

  # Create path?
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





###########################################################################/**
# @RdocMethod updateSamplesFile
#
# @title "Updates the samples.js file"
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
#  Returns (invisibly) the pathname to the samples.js file.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("updateSamplesFile", "ChromosomeExplorer", function(this, ..., verbose=FALSE) {
  require("R.rsp") || throw("Package not loaded: R.rsp");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  path <- getPath(this);
  parentPath <- getParent(path);
  parent2Path <- getParent(parentPath);
  parent3Path <- getParent(parent2Path);

  verbose && enter(verbose, "Compiling samples.js");
  srcPath <- getTemplatePath(this);
  pathname <- filePath(srcPath, "rsp", "ChromosomeExplorer", "samples.js.rsp");
  verbose && cat(verbose, "Source: ", pathname);
  outFile <- gsub("[.]rsp$", "", basename(pathname));
  outPath <- parent3Path;
  verbose && cat(verbose, "Output path: ", outPath);


  verbose && enter(verbose, "Scanning directories for available chip types");
  # Find all directories matching these
  dirs <- list.files(path=parent2Path, full.names=TRUE);
  dirs <- dirs[sapply(dirs, FUN=isDirectory)];

  # Get possible chip types
  model <- getModel(this);
  chipTypes <- c(getChipType(model), getChipTypes(model));
  chipTypes <- intersect(chipTypes, basename(dirs));
  verbose && cat(verbose, "Detected chip types: ", 
                                           paste(chipTypes, collapse=", "));
  verbose && exit(verbose);

  # Get available zooms
  verbose && enter(verbose, "Scanning image files for available zooms");
  pattern <- ".*,x([0-9][0-9]*)[.]png$";
  zooms <- list.files(path=path, pattern=pattern);
  zooms <- gsub(pattern, "\\1", zooms);
  zooms <- gsub("^0*", "", zooms);
  if (length(zooms) == 0) {
    # Default zooms
    zooms <- c(1, 2, 4, 8, 16, 32, 64, 128);
  }
  zooms <- unique(zooms);
  zooms <- as.integer(zooms);
  zooms <- sort(zooms);
  verbose && cat(verbose, "Detected (or default) zooms: ", paste(zooms, collapse=", "));
  verbose && exit(verbose);

  # Compile RSP file
  verbose && enter(verbose, "Compiling RSP");
  env <- new.env();
  env$chipTypes <- chipTypes;
  env$samples <- getFullNames(this);
  env$sampleLabels <- getNames(this);
  env$zooms <- zooms;
  pathname <- rspToHtml(pathname, path=NULL, 
                        outFile=outFile, outPath=outPath, 
                        overwrite=TRUE, envir=env);
  verbose && exit(verbose);


  verbose && exit(verbose);
  
  invisible(pathname);
}, private=TRUE)



setMethodS3("addIndexFile", "ChromosomeExplorer", function(this, filename="index.html", ...) {
  NextMethod("addIndexFile", this, filename=filename, ...);
}, protected=TRUE)




setMethodS3("setup", "ChromosomeExplorer", function(this, ..., force=FALSE) {
  # Setup includes/?
  addIncludes(this, ..., force=force);

  # Setup index.html
  addIndexFile(this, ..., force=force);

  # Setup samples.js?
  updateSamplesFile(this, ...);
}, private=TRUE)



setMethodS3("writeGraphs", "ChromosomeExplorer", function(x, arrays=NULL, ...) {
  # To please R CMD check.
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  if (is.null(arrays))
    arrays <- getNames(this);


  # Get the model
  model <- getModel(this);

  path <- getPath(this);
  plotband <- this$.plotCytoband;  # Plot cytoband?
  plot(model, path=path, imageFormat="png", plotband=plotband, arrays=arrays, ...);

  invisible(path);
}, private=TRUE)



setMethodS3("writeRegions", "ChromosomeExplorer", function(this, arrays=NULL, nbrOfSnps=c(3,Inf), smoothing=c(-Inf,-0.15, +0.15,+Inf), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  if (is.null(arrays)) 
    arrays <- getNames(this);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Writing CN regions");

  # Extract and write regions
  model <- getModel(this);
  pathname <- writeRegions(model, arrays=arrays, nbrOfSnps=nbrOfSnps, smoothing=smoothing, ..., skip=FALSE, verbose=less(verbose));

  dest <- filePath(getPath(this), "regions.xls");
  res <- file.copy(pathname, dest, overwrite=TRUE);
  if (!res)
    dest <- NULL;

  verbose && exit(verbose);

  invisible(dest);
}, private=TRUE)



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
#   \item{chromosome}{A @vector of chromosomes specifying which
#     chromosomes to be considered.  If @NULL, all are processed.}
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
setMethodS3("process", "ChromosomeExplorer", function(this, arrays=NULL, chromosomes=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  if (is.null(arrays))
    arrays <- getNames(this);

  # Argument 'chromosomes':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  
  verbose && enter(verbose, "Generating ChromosomeExplorer report");

  # Setup HTML, CSS, Javascript files first
  setup(this, ..., verbose=less(verbose));

  # Generate bitmap images
  writeGraphs(this, arrays=arrays, chromosomes=chromosomes, ..., verbose=less(verbose));

  # Update samples.js
  updateSamplesFile(this, ..., verbose=less(verbose));

  # Write regions file
  writeRegions(this, arrays=arrays, chromosomes=chromosomes, ..., verbose=less(verbose));

  verbose && exit(verbose);
})


setMethodS3("display", "ChromosomeExplorer", function(this, filename="index.html", ...) {
  NextMethod("display", this, filename=filename, ...);
})


##############################################################################
# HISTORY:
# 2007-05-08
# o Added default zoom levels to updateSamplesFile() for ChromosomeExplorer.  
#   This is applies the first time process() is called.
# 2007-03-19
# o Now ChromosomeExplorer extends Explorer.
# 2007-02-06
# o Now templates are in reports/templates/ and includes in reports/includes/.
# o Updated the path to <rootPath>/<dataSetName>/<tags>/<chipType>/<set>/.
# 2007-01-17
# o Now all 'arrays' arguments can contain array names.
# o Added getArrays() and setArrays() in order to focus on a subset of the
#   arrays in the model.
# 2007-01-16
# o Added getAlias() and setAlias(), and updated getName() accordingly.
#   This makes it easy to change the name of output set for subsets of
#   arrays.
# 2007-01-15
# o Added some more Rdoc comments.
# 2007-01-10
# o BUG FIX: setup() would only add index.html if includes/ were missing.
# 2007-01-08
# o Added display().
# 2007-01-07
# o TODO: The region filter for writeRegions() is hardwired for now.
# o Update to include regions.xls file on the CE output page too.
# o Now all HTML, CSS, and Javascript files are created too.
# 2007-01-04
# o Created.
##############################################################################
