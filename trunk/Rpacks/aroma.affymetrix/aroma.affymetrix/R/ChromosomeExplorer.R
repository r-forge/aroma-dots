###########################################################################/**
# @RdocClass ChromosomeExplorer
#
# @title "The ChromosomeExplorer class"
#
# \description{
#  @classhierarchy
#
#  This class represents a Chromosome Explorer report.
#  It provides methods to define the report and generate it.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("ChromosomeExplorer", function(model=NULL, arrays=NULL, chromosomes=NULL, ...) {
  this <- extend(Object(), "ChromosomeExplorer",
    .model = NULL,
    .arrays = NULL,
    .chromosomes = NULL
  );

  if (!is.null(model)) {
    setModel(this, model);
    setArrays(this, arrays);
    setChromosomes(this, chromosomes);
  }

  this;
})


setMethodS3("getModel", "ChromosomeExplorer", function(this, ...) {
  this$.model;
})

setMethodS3("setModel", "ChromosomeExplorer", function(this, model, ...) {
  # Arguments 'model':
  if (!inherits(model, "GladModel")) {
    throw("Argument 'model' must be a GladModel: ", class(model)[1]);
  }

  this$.model <- model;
})

setMethodS3("getArrays", "ChromosomeExplorer", function(this, ...) {
  arrays <- this$.arrays;
  if (is.null(arrays))
    arrays <- seq(length=nbrOfArrays(getModel(this)));
  arrays;
})

setMethodS3("setArrays", "ChromosomeExplorer", function(this, arrays=NULL, ...) {
  # Arguments 'arrays':
  if (!is.null(arrays)) {
    arrays <- Arguments$getIndices(arrays, 
                                  range=c(1, nbrOfArrays(getModel(this))));
  }

  this$.arrays <- arrays;
})


setMethodS3("getChromosomes", "ChromosomeExplorer", function(this, ...) {
  chromosomes <- this$.chromosomes;
  if (is.null(chromosomes))
    chromosomes <- 1:23;
  chromosomes;
})

setMethodS3("setChromosomes", "ChromosomeExplorer", function(this, chromosomes=NULL, ...) {
  # Arguments 'chromosomes':
  if (!is.null(chromosomes)) {
    if (is.character(chromosomes))
      chromosomes <- match(chromosomes, c(1:22,"X","Y")); 
    chromosomes <- Arguments$getIndices(chromosomes, range=c(1,23)); 
  }

  this$.chromosomes <- chromosomes;
})


setMethodS3("getPath", "ChromosomeExplorer", function(this, ...) {
  rootPath <- "chromosomeExplorer";
  model <- getModel(this);
  name <- getName(model);
  path <- filePath(rootPath, name);
  mkdirs(path);
  path;
})


setMethodS3("getFigurePath", "ChromosomeExplorer", function(this, ...) {
  path <- getPath(this, ...);
  model <- getModel(this);
  name <- getName(model);
  tags <- getTags(model);
  cdf <- getCdf(model);
  chipType <- getChipType(cdf);
  chipType <- gsub("-monocell", "", chipType);
  tags <- paste(tags, collapse=",");
  path <- filePath(path, tags, chipType, expandLinks="any");
  mkdirs(path);
  path;
})


setMethodS3("generatePlots", "ChromosomeExplorer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate all graphs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  model <- getModel(this);
  arrays <- getArrays(this);
  chromosomes <- getChromosomes(this);

  path <- getFigurePath(this);

  verbose && enter(verbose, "Generating all plots");
  verbose && cat(verbose, "Arrays: ", seqToHumanReadable(arrays));
  verbose && cat(verbose, "Chromosomes: ", seqToHumanReadable(chromosomes));
  verbose && cat(verbose, "Figure path: ", path);

  plot(model, arrays=arrays, chromosomes=chromosomes, path=path, ..., 
                                                       verbose=less(verbose));
  verbose && exit(verbose);
}, protected=TRUE)


setMethodS3("generate", "ChromosomeExplorer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Generating ChromosomeExplorer report");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  model <- getModel(this);
  verbose && cat(verbose, "Model:");
  verbose && print(verbose, model);

  # Create paths

  verbose && enter(verbose, "Copying files and directories");
  srcPath <- "~/braju.com.R/aroma.affymetrix/aroma.affymetrix/inst/chromosomeExplorer/";
  path <- getPath(this);
  for (file in c("index.html")) {
    from <- filePath(srcPath, file);
    to <- filePath(path, file);
    file.copy(from, to);
  }

  for (dir in c("css", "images", "js")) {
    from <- filePath(srcPath, dir);
    to <- filePath(path, dir);
    copyDirectory(from=from, to=to);
  }
  verbose && exit(verbose);

  # Generate plots
#  generatePlots(this, ..., verbose=verbose);


  verbose && exit(verbose);
})


############################################################################
# HISTORY:
# 2006-11-24
# o Recreated again.
# 2006-09-11
# o Recreated.
############################################################################
