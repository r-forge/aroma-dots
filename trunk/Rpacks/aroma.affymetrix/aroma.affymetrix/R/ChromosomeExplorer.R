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
  env$samples <- getArrays(model);
  env$zooms <- zooms;
  pathname <- rspToHtml(pathname, path=NULL, 
                        outFile=outFile, outPath=outPath, 
                        overwrite=TRUE, envir=env);
  verbose && exit(verbose);


  verbose && exit(verbose);
  
  invisible(pathname);
}, private=TRUE)


setMethodS3("addIncludes", "ChromosomeExplorer", function(this, ..., overwrite=FALSE, verbose=FALSE) {
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
  if (!isDirectory(destPath)) {
    verbose && enter(verbose, "Copying template files");
    verbose && cat(verbose, "Source path: ", srcPath);
    verbose && cat(verbose, "Destination path: ", destPath);
    pathnames <- copyDirectory(from=srcPath, to=destPath, recursive=TRUE);
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Copying index.html");
  srcPathname <- filePath(srcPath, "html", "index.html");
  outPath <- getParent(getPath(this));
  outPathname <- filePath(outPath, "index.html");
  verbose && cat(verbose, "Source pathname: ", srcPathname);
  verbose && cat(verbose, "Destination pathname: ", outPathname);
  file.copy(srcPathname, outPathname, overwrite=overwrite);
  verbose && exit(verbose);

  verbose && exit(verbose);
}, private=TRUE)


setMethodS3("setup", "ChromosomeExplorer", function(this, ..., force=FALSE) {
  # Setup includes/?
  path <- filePath(getRootPath(this), "includes");
  if (force || !isDirectory(path)) {
    addIncludes(this, ...);
  }

  # Setup samples.js?
  updateSamplesFile(this, ...);
}, private=TRUE)


setMethodS3("writeGraphs", "ChromosomeExplorer", function(x, ...) {
  # To please R CMD check.
  this <- x;

  path <- getPath(this);
  model <- getModel(this);
  plot(model, path=path, imageFormat="png", ...);

  invisible(path);
}, private=TRUE)


setMethodS3("writeRegions", "ChromosomeExplorer", function(this, nbrOfSnps=c(3,Inf), smoothing=c(-Inf,-0.15, +0.15,+Inf), ...) {
  # Extract and write regions
  model <- getModel(this);
  pathname <- writeRegions(model, nbrOfSnps=nbrOfSnps, smoothing=smoothing, ..., skip=FALSE);

  dest <- filePath(getPath(this), "regions.xls");
  res <- file.copy(pathname, dest, overwrite=TRUE);
  if (!res)
    dest <- NULL;
  invisible(dest);
}, private=TRUE)



setMethodS3("process", "ChromosomeExplorer", function(this, arrays=NULL, chromosomes=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':

  # Argument 'chromosomes':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Generating ChromosomeExplorer report");

  # Setup HTML, CSS, Javascript files first
  setup(this, ...);

  # First the model, just in case
  model <- getModel(this);
  fit(model, arrays=arrays, chromosomes=chromosomes, ...);

  # Generate bitmap images
  writeGraphs(this, arrays=arrays, chromosomes=chromosomes, ...);

  # Update samples.js
  updateSamplesFile(this, ...);

  # Write regions file
  writeRegions(this, arrays=arrays, chromosomes=chromosomes, ...);

  verbose && exit(verbose);
})


##############################################################################
# HISTORY:
# 2007-01-07
# o TODO: The region filter for writeRegions() is hardwired for now.
# o Update to include regions.xls file on the CE output page too.
# o Now all HTML, CSS, and Javascript files are created too.
# 2007-01-04
# o Created.
##############################################################################
