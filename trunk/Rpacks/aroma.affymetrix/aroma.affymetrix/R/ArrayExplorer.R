###########################################################################/**
# @RdocClass ArrayExplorer
#
# @title "The ArrayExplorer class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{csTuple}{An @see "AffymetrixCelSet" object.}
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
setConstructorS3("ArrayExplorer", function(csTuple=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csTuple':
  if (!is.null(csTuple)) {
    if (!inherits(csTuple, "AffymetrixCelSetTuple")) {
      csTuple <- AffymetrixCelSetTuple(csTuple);
    }
  }

  extend(Explorer(...), "ArrayExplorer",
    .csTuple = csTuple,
    .reporters = NULL
  )
})



setMethodS3("as.character", "ArrayExplorer", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Number of chip types: %d", nbrOfChipTypes(this)));
  s <- c(s, paste("Number of arrays:", nbrOfArrays(this)));
  colorMaps <- getColorMaps(this);
  if (length(colorMaps) == 0) {
    colorMaps <- "<not specified>";
  } else {
    colorMaps <- paste(colorMaps, collapse="; ")
  }
  s <- c(s, paste("Color maps:", colorMaps));
  s <- c(s, sprintf("Main path: %s", getMainPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getListOfReporters", "ArrayExplorer", function(this, ...) {
  reporters <- this$.reporters;

  if (is.null(reporters)) {
    tags <- getTags(this);
    setTuple <- getSetTuple(this);
    csList <- getListOfSets(setTuple);
    reporters <- lapply(csList, FUN=function(cs) {
      reporter <- SpatialReporter(cs, tags=tags);
      reporter;
    });
    this$.reporters <- reporters;
  }

  reporters;
}, protected=TRUE);


setMethodS3("setAlias", "ArrayExplorer", function(this, ...) {
  NextMethod("setAlias", this, ...);
  reporters <- getListOfReporters(this);
  lapply(reporters, FUN=setAlias, ...);
  invisible(this);
})

setMethodS3("getAlias", "ArrayExplorer", function(this, ...) {
  reporters <- getListOfReporters(this);
  getAlias(reporters[[1]], ...);
})


setMethodS3("getAsteriskTags", "ArrayExplorer", function(this, ...) {
  "";
})



###########################################################################/**
# @RdocMethod getDataSet
#
# @title "Gets the data set"
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
#  Returns a @see "AffymetrixCelSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getDataSet", "ArrayExplorer", function(this, ...) {
  csTuple <- getSetTuple(this);
  csList <- getListOfSets(csTuple);
  csList[[1]];  # AD HOC for now. /HB 2007-03-19
})

setMethodS3("getSetTuple", "ArrayExplorer", function(this, ...) {
  this$.csTuple;
})

setMethodS3("getNameOfInput", "ArrayExplorer", function(this, ...) {
  st <- getSetTuple(this);
  getName(st, ...);
}, protected=TRUE)


setMethodS3("getTagsOfInput", "ArrayExplorer", function(this, ...) {
  st <- getSetTuple(this);
  getTags(st, ...);
}, protected=TRUE)


setMethodS3("nbrOfChipTypes", "ArrayExplorer", function(this, ...) {
  nbrOfChipTypes(getSetTuple(this));
})

setMethodS3("getListOfCdfs", "ArrayExplorer", function(this, ...) {
  getListOfCdfs(getSetTuple(this));
}, private=TRUE)


setMethodS3("getArraysOfInput", "ArrayExplorer", function(this, ...) {
  setTuple <- getSetTuple(this);
  getArrays(setTuple);
}, protected=TRUE)



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
setMethodS3("setArrays", "ArrayExplorer", function(this, arrays=NULL, ...) {
  # Argument 'arrays':
  if (!is.null(arrays)) {
    setTuple <- getSetTuple(this);
    arrays <- indexOfArrays(setTuple, arrays);
    arrays <- getArrays(setTuple)[arrays];
  }

  this$.arrays <- arrays;
})



setMethodS3("updateOnChipTypeJS", "ArrayExplorer", function(this, ..., verbose=FALSE) {
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

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mainPath <- getMainPath(this);
  setTuple <- getSetTuple(this);

  outFile <- sprintf("%s.onChipType.js", class(this)[1]);
  filename <- sprintf("%s.rsp", outFile);
  verbose && enter(verbose, "Updating ", outFile);
  srcPath <- getTemplatePath(this);
  pathname <- filePath(srcPath, "rsp", class(this)[1], filename);
  verbose && cat(verbose, "Source: ", pathname);

  outPath <- mainPath;
  verbose && cat(verbose, "Output path: ", outPath);

  # For each SpatialReporters
  reporters <- getListOfReporters(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Compile RSP file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Compiling RSP");
  env <- new.env();
  env$reporters <- reporters;
  pathname <- rspToHtml(pathname, path=NULL, 
                        outFile=outFile, outPath=outPath, 
                        overwrite=TRUE, envir=env);
  verbose && exit(verbose);


  verbose && exit(verbose);
  
  invisible(pathname);
}, private=TRUE)





setMethodS3("updateOnLoadJS", "ArrayExplorer", function(this, ..., verbose=FALSE) {
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mainPath <- getMainPath(this);
  setTuple <- getSetTuple(this);

  outFile <- sprintf("%s.onLoad.js", class(this)[1]);
  filename <- sprintf("%s.rsp", outFile);

  verbose && enter(verbose, "Updating ", outFile);
  srcPath <- getTemplatePath(this);
  pathname <- filePath(srcPath, "rsp", class(this)[1], filename);
  verbose && cat(verbose, "Source: ", pathname);

  outPath <- mainPath;
  verbose && cat(verbose, "Output path: ", outPath);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Compile RSP file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipTypes <- getChipTypes(setTuple);
  verbose && cat(verbose, "Detected chip types: ", 
                                           paste(chipTypes, collapse=", "));
  verbose && enter(verbose, "Compiling RSP");
  env <- new.env();
  env$chipTypes <- chipTypes;
  arrays <- getFullNames(setTuple);
  pathname <- rspToHtml(pathname, path=NULL, 
                        outFile=outFile, outPath=outPath, 
                        overwrite=TRUE, envir=env);
  verbose && exit(verbose);

  verbose && exit(verbose);
  
  invisible(pathname);
}, private=TRUE)



setMethodS3("setup", "ArrayExplorer", function(this, ..., force=FALSE) {
  # Setup includes/
  addIncludes(this, ..., force=force);

  # Setup HTML explorer page
  addIndexFile(this, ..., force=force);

  # Update Javascript files
  updateOnChipTypeJS(this, ...);
  updateOnLoadJS(this, ...);
}, private=TRUE)


setMethodS3("addColorMap", "ArrayExplorer", function(this, ...) {
  reporters <- getListOfReporters(this);
  res <- lapply(reporters, FUN=addColorMap, ...);
  invisible(res);
})


setMethodS3("setColorMaps", "ArrayExplorer", function(this, ...) {
  reporters <- getListOfReporters(this);
  res <- lapply(reporters, FUN=setColorMaps, ...);
  invisible(res);
})

setMethodS3("getColorMaps", "ArrayExplorer", function(this, parsed=FALSE, ...) {
  # All reporters should have the same color maps
  reporter <- getListOfReporters(this)[[1]];
  getColorMaps(reporter, ...);
})


setMethodS3("setScale", "ArrayExplorer", function(this, ...) {
  reporters <- getListOfReporters(this);
  lapply(reporters, FUN=setScale, ...);
  invisible(this);
})

setMethodS3("getScale", "ArrayExplorer", function(this, parsed=FALSE, ...) {
  # All reporters should have the same color maps
  reporter <- getListOfReporters(this)[[1]];
  getScale(reporter, ...);
})


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
setMethodS3("process", "ArrayExplorer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  
  verbose && enter(verbose, "Generating ", class(this)[1], " report");

  # Setup HTML, CSS, Javascript files first
  setup(this, ..., verbose=less(verbose));

  # Generate bitmap images
  reporters <- getListOfReporters(this);
  res <- lapply(reporters, FUN=function(reporter) {
    writeImages(reporter, aliases=names(getArrays(this)), ..., verbose=less(verbose));
  });

  # Update Javascript files
  updateOnChipTypeJS(this, ..., verbose=less(verbose));
  updateOnLoadJS(this, ..., verbose=less(verbose));

  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);
})




##############################################################################
# HISTORY:
<<<<<<< .mine
# 2007-03-30
# o Added getScale() and setScale().
=======
# 2007-08-09
# o Renamed updateSampleFile() to updateOnLoadJS().
# o Added updateOnChipTypeJS().
>>>>>>> .r2966
# 2007-03-20
# o Removed argument arrays from process().
# o Added setAlias() which also sets the alias on the reporters.
# o Added getAlias() to inherit the alias from the reporters.
# 2007-03-19
# o Class can now handle multiple chip types.
# o Class is now making use of the SpatialReporter class.
# o Now ChromosomeExplorer extends Explorer.
# o BUG FIX: setArrays() called indexOfArrays() instead of indexOf().
# 2007-02-28
# o BUG FIX: setColorMaps() gave "Error in addColorMap.ArrayExplorer(this, 
#   colorMap, ...) : object "nbrOfColors" not found".
# 2007-02-08
# o Created from ChromosomeExplorer.R.
##############################################################################
