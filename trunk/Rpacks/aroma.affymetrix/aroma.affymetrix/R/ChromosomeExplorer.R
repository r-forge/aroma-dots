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
    .alias = NULL,
    .model = model,
    .arrays = NULL,
    .plotCytoband = TRUE
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
setMethodS3("getArrays", "ChromosomeExplorer", function(this, ...) {
  arrays <- this$.arrays;
  if (is.null(arrays)) {
    model <- getModel(this);
    arrays <- getArrays(model);
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
    model <- getModel(this);
    arrays <- indexOfArrays(model, arrays=arrays);
    arrays <- getArrays(model)[arrays];
  }

  this$.arrays <- arrays;
})





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
setMethodS3("nbrOfArrays", "ChromosomeExplorer", function(this, ...) {
  arrays <- getArrays(this);
  length(arrays);
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
setMethodS3("getAlias", "ChromosomeExplorer", function(this, ...) {
  this$.alias;
})



###########################################################################/**
# @RdocMethod setAlias
#
# @title "Sets the alias of the output set"
#
# \description{
#   @get "title".
#   If specified, the alias overrides the model name, which is used by 
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
setMethodS3("setAlias", "ChromosomeExplorer", function(this, alias=NULL, ...) {
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
#  @get "title", which is the same as the name of the model.
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
#  If a name alias has not been set explicitly, the name of the model will
#  used.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "ChromosomeExplorer", function(this, ...) {
  name <- getAlias(this);
  if (is.null(name)) {
    model <- getModel(this);
    name <- getName(model, ...);
  }
  name;
})



###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the explorer"
#
# \description{
#  @get "title", which are the tags of the model plus additional tags.
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


setMethodS3("getTemplatePath", "ChromosomeExplorer", function(this, ..., verbose=FALSE) {
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
  path <- filePath("ce/.template", expandLinks="any");
  if (!isDirectory(path)) {
    path <- system.file("chromosomeExplorer", "includes", 
                                            package="aroma.affymetrix");
  }
  verbose && exit(verbose);

  path;
}, private=TRUE)




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
  require(R.rsp) || throw("Package not loaded: R.rsp");

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

  verbose && enter(verbose, "Compiling samples.js");
  srcPath <- getTemplatePath(this);
  pathname <- filePath(srcPath, "rsp", "samples.js.rsp");
  verbose && cat(verbose, "Source: ", pathname);
  outFile <- gsub("[.]rsp$", "", basename(pathname));
  outPath <- getParent(parentPath);
  verbose && cat(verbose, "Output path: ", outPath);


  verbose && enter(verbose, "Scanning directories for available chip types");
  # Find all directories matching these
  dirs <- list.files(path=parentPath, full.names=TRUE);
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
  zooms <- unique(zooms);
  zooms <- as.integer(zooms);
  zooms <- sort(zooms);
  verbose && cat(verbose, "Detected zooms: ", paste(zooms, collapse=", "));
  verbose && exit(verbose);

  # Compile RSP file
  verbose && enter(verbose, "Compiling RSP");
  env <- new.env();
  env$chipTypes <- chipTypes;
  env$samples <- getArrays(this);
  env$zooms <- zooms;
  pathname <- rspToHtml(pathname, path=NULL, 
                        outFile=outFile, outPath=outPath, 
                        overwrite=TRUE, envir=env);
  verbose && exit(verbose);


  verbose && exit(verbose);
  
  invisible(pathname);
}, private=TRUE)




setMethodS3("addIncludes", "ChromosomeExplorer", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Setting up ChromsomeExplorer report files");

  srcPath <- getTemplatePath(this);

  destPath <- filePath(getRootPath(this), "includes");
  if (force || !isDirectory(destPath)) {
    verbose && enter(verbose, "Copying template files");
    verbose && cat(verbose, "Source path: ", srcPath);
    verbose && cat(verbose, "Destination path: ", destPath);
    pathnames <- copyDirectory(from=srcPath, to=destPath, recursive=TRUE);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);
})

setMethodS3("addIndexFile", "ChromosomeExplorer", function(this, ..., force=FALSE, verbose=FALSE) {
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
  srcPathname <- filePath(srcPath, "html", "index.html");
  outPath <- getParent(getPath(this));
  outPathname <- filePath(outPath, "index.html");

  if (force || !isFile(outPathname)) {
    verbose && enter(verbose, "Copying index.html");
    verbose && cat(verbose, "Source pathname: ", srcPathname);
    verbose && cat(verbose, "Destination pathname: ", outPathname);
    file.copy(srcPathname, outPathname, overwrite=TRUE);
    verbose && exit(verbose);
  }
}, private=TRUE)




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
    arrays <- getArrays(this);



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
    arrays <- getArrays(this);

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
    arrays <- getArrays(this);

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
setMethodS3("display", "ChromosomeExplorer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Opening ChromosomeExplorer");

  # The path to the explorer HTML document
  path <- getPath(this);
  path <- getParent(path);
  pathname <- filePath(path, "index.html", expandLinks="any");

  # Just in case, is setup needed?
  if (!isFile(pathname)) {
    setup(this, verbose=less(verbose));
    if (!isFile(pathname)) {
      throw("Cannot open ChromosomeExplorer. No such file: ", pathname);
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
