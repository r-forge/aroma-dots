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
#   \item{celSet}{An @see "AffymetrixCelSet" object.}
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
setConstructorS3("ArrayExplorer", function(celSet=NULL, tags="*", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'celSet':
  if (!is.null(celSet)) {
    if (!inherits(celSet, "AffymetrixCelSet")) {
      throw("Argument 'celSet' is not a 'AffymetrixCelSet': ", class(celSet)[1]);
    }
  }

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    tags[tags == "*"] <- "";
  }


  extend(Object(), "ArrayExplorer",
    .alias = NULL,
    .celSet = celSet,
    .arrays = NULL
  )
})


setMethodS3("as.character", "ArrayExplorer", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, paste("Number of arrays:", nbrOfArrays(this)));
  s <- c(s, paste("Color maps:", paste(getColorMaps(this), collapse="; ")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)



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
  this$.celSet;
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
setMethodS3("getArrays", "ArrayExplorer", function(this, ...) {
  arrays <- this$.arrays;
  if (is.null(arrays)) {
    celSet <- getDataSet(this);
    arrays <- sapply(celSet, getFullName);
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
setMethodS3("setArrays", "ArrayExplorer", function(this, arrays=NULL, ...) {

  # Argument 'arrays':
  if (!is.null(arrays)) {
    celSet <- getDataSet(this);
    arrays <- indexOfArrays(celSet, arrays=arrays);
    arrays <- sapply(celSet, getFullName)[arrays];
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
setMethodS3("nbrOfArrays", "ArrayExplorer", function(this, ...) {
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
setMethodS3("getAlias", "ArrayExplorer", function(this, ...) {
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
setMethodS3("setAlias", "ArrayExplorer", function(this, alias=NULL, ...) {
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
setMethodS3("getName", "ArrayExplorer", function(this, ...) {
  name <- getAlias(this);
  if (is.null(name)) {
    celSet <- getDataSet(this);
    name <- getName(celSet, ...);
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
setMethodS3("getTags", "ArrayExplorer", function(this, ...) {
  celSet <- getDataSet(this);
  tags <- getTags(celSet);
  tags <- c(tags, this$.tags);
  tags <- unique(tags);
  tags;
})


setMethodS3("getFullName", "ArrayExplorer", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "ArrayExplorer", function(this, ...) {
  "reports";
}, private=TRUE)


setMethodS3("getPath", "ArrayExplorer", function(this, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  name <- getName(this);

  # Tags
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  if (nchar(tags) == 0) {
    tags <- "raw";  # Default
  }

  # Chip type    
  celSet <- getDataSet(this);
  cdf <- getCdf(celSet);
  chipType <- getChipType(cdf);

  # Image set
  set <- "spatial";

  # The full path
  path <- filePath(rootPath, name, tags, chipType, set, expandLinks="any");
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
})



setMethodS3("getTemplatePath", "ArrayExplorer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating template files for ", class(this)[1]);
  # Search for template files
  rootPath <- getRootPath(this);
  path <- filePath(rootPath, "templates", expandLinks="any");
  if (!isDirectory(path)) {
    path <- system.file("reports", "templates", package="aroma.affymetrix");
  }
  verbose && exit(verbose);

  path;
}, private=TRUE)


setMethodS3("getIncludePath", "ArrayExplorer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating include files for ", class(this)[1]);
  # Search for include files
  path <- system.file("reports", "includes", package="aroma.affymetrix");
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
setMethodS3("updateSamplesFile", "ArrayExplorer", function(this, ..., verbose=FALSE) {
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

  arrays <- getArrays(this);

  path <- getPath(this);
  parentPath <- getParent(path);
  parent2Path <- getParent(parentPath);

  outFile <- sprintf("%s.onLoad.js", class(this)[1]);
  filename <- sprintf("%s.rsp", outFile);
  verbose && enter(verbose, "Compiling ", outFile);
  srcPath <- getTemplatePath(this);
  pathname <- filePath(srcPath, "rsp", class(this)[1], filename);
  verbose && cat(verbose, "Source: ", pathname);
  outPath <- getParent(parent2Path);
  verbose && cat(verbose, "Output path: ", outPath);


  verbose && enter(verbose, "Scanning directories for available chip types");
  # Find all directories matching these
  dirs <- list.files(path=parent2Path, full.names=TRUE);
  dirs <- dirs[sapply(dirs, FUN=isDirectory)];
  # Get possible chip types
  celSet <- getDataSet(this);
  cdf <- getCdf(celSet);
  chipTypes <- c(getChipType(cdf));
  chipTypes <- intersect(chipTypes, basename(dirs));
  verbose && cat(verbose, "Detected chip types: ", 
                                           paste(chipTypes, collapse=", "));
  verbose && exit(verbose);

  # Get available color maps
  verbose && enter(verbose, "Scanning image files for available color maps");
  pattern <- "[^,]*,.*[.]png$";
  colorMaps <- list.files(path=path, pattern=pattern);
  colorMaps <- gsub("[.]png$", "", colorMaps);
  for (array in arrays) {
    pattern <- sprintf("^%s,", array);
    colorMaps <- gsub(pattern, "", colorMaps);
  }
  colorMaps <- unique(colorMaps);
  colorMaps <- sort(colorMaps);
  verbose && cat(verbose, "Detected color maps: ", paste(colorMaps, collapse=", "));
  verbose && exit(verbose);

  # Compile RSP file
  verbose && enter(verbose, "Compiling RSP");
  env <- new.env();
  env$chipTypes <- chipTypes;
  env$aliases <- gsub(",.*", "", arrays);
  env$samples <- arrays; 
  env$colorMaps <- colorMaps;
  pathname <- rspToHtml(pathname, path=NULL, 
                        outFile=outFile, outPath=outPath, 
                        overwrite=TRUE, envir=env);
  verbose && exit(verbose);


  verbose && exit(verbose);
  
  invisible(pathname);
}, private=TRUE)




setMethodS3("addIncludes", "ArrayExplorer", function(this, ..., force=FALSE, verbose=FALSE) {
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

setMethodS3("addIndexFile", "ArrayExplorer", function(this, ..., force=FALSE, verbose=FALSE) {
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
  filename <- sprintf("%s.html", class(this)[1]);
  srcPathname <- filePath(srcPath, "html", class(this)[1], filename);
  outPath <- getParent(getParent(getPath(this)));
  outPathname <- filePath(outPath, filename);

  if (force || !isFile(outPathname)) {
    verbose && enter(verbose, "Copying ", filename);
    verbose && cat(verbose, "Source pathname: ", srcPathname);
    verbose && cat(verbose, "Destination pathname: ", outPathname);
    if (!isFile(srcPathname))
      throw("File not found: ", srcPathname);
    file.copy(srcPathname, outPathname, overwrite=TRUE);
    verbose && exit(verbose);
  }
}, private=TRUE)




setMethodS3("setup", "ArrayExplorer", function(this, ..., force=FALSE) {
  # Setup includes/?
  addIncludes(this, ..., force=force);

  # Setup index.html
  addIndexFile(this, ..., force=force);

  # Setup samples.js?
  updateSamplesFile(this, ...);
}, private=TRUE)


setMethodS3("addColorMap", "ArrayExplorer", function(this, colorMap, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'colorMap':
  colorMap <- Arguments$getCharacter(colorMap, nchar=c(1,Inf));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Parser argument 'colorMap'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  tags <- colorMap;
  parts <- strsplit(tags, split=",")[[1]];
  n <- length(parts);
  if (n > 2) {
    throw("Argument 'colorMap' must not contain more than two tags: ", tags);
  }

  if (n == 1) {
    transform <- NULL;
    palette <- parts[1];
  } else if (n == 2) {
    transform <- parts[1];
    palette <- parts[2];
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Setup color map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.null(transform)) {
    transform <- sqrt;
  } else if (is.character(transform)) {
    # Check transform
    if (!exists(transform, mode="function")) {
      throw("Argument 'colorMap' specifies an unknown transform function ('",
                                               transform, "'): ", colorMap);
    }
    transform <- get(transform, mode="function");
  } else if (!is.function(transform)) {
    throw("Argument 'colorMap' specifies an invalid transform: ", 
                                                           mode(transform));
  }

  if (is.null(nbrOfColors)) {
    nbrOfColors <- 256;
  } else if (is.na(nbrOfColors)) {
    throw("Argument 'colorMap' specifies an invalid number of colors ('", 
                                                   ncol, "'): ", colorMap);
  }

  if (is.null(palette)) {
    palette <- gray.colors(256);
  } else if (is.function(palette)) {
    palette <- palette(256);
  } else if (is.character(palette)) {
    # Parse color palette tag
    pattern <- "^([^_]*)(|[_][0-9]*)$";
    name <- gsub(pattern, "\\1", palette);
    nbrOfColors <- gsub(pattern, "\\2", palette);
    nbrOfColors <- gsub("[_]", "", nbrOfColors);
    nbrOfColors <- as.integer(nbrOfColors);
    if (is.na(nbrOfColors))
      nbrOfColors <- 256;

    # Search for function <palette>() and then <palette>.colors()
    if (!exists(name, mode="function")) {
      name <- sprintf("%s.colors", palette);
      if (!exists(name, mode="function")) {
        throw("Argument 'colorMap' specifies an unknown palette function ('", 
                                                   name, "'): ", colorMap);
      }
    }
    fcn <- get(name, mode="function");
    palette <- fcn(nbrOfColors);
  }
  

  map <- list(list(
    tags = tags,
    transform = transform,
    palette = palette
  ));
  names(map) <- map$tags;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Add color map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  colorMaps <- this$.colorMaps;
  if (is.null(colorMaps))
    colorMaps <- list();  
  colorMaps <- c(colorMaps, map);
  this$.colorMaps <- colorMaps;
})


setMethodS3("setColorMaps", "ArrayExplorer", function(this, colorMaps=c("sqrt,yellow", "sqrt,rainbow"), ...) {
  this$.colorMaps <- NULL;
  for (colorMap in colorMaps) {
    addColorMap(this, colorMap, ...);
  }
})

setMethodS3("getColorMaps", "ArrayExplorer", function(this, parsed=FALSE, ...) {
  colorMaps <- this$.colorMaps;

  if (!parsed)
    colorMaps <- sapply(colorMaps, .subset2, "tags");

  colorMaps;
})

setMethodS3("writeImages", "ArrayExplorer", function(this, arrays=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  if (is.null(arrays))
    arrays <- getArrays(this);



  # Get the path to the image directory
  path <- getPath(this);

  # Get the CEL set of interest
  cs <- getDataSet(this);
  if (!is.null(arrays)) {
    allArrays <- sapply(cs, getFullName);
    idxs <- match(arrays, allArrays);
    cs <- extract(cs, idxs);
  }

  # Get the color maps to be generated
  colorMaps <- getColorMaps(this, parsed=TRUE);
  if (length(colorMaps) == 0) {
    warning("No color maps specified. Nothing to do.");
    return(invisible(path));
  }

  # For each array...
  for (kk in seq(cs)) {
    df <- getFile(cs, kk);
    verbose && enter(verbose, sprintf("Array #%d ('%s')", kk, getName(df)));
    # For each color map...
    for (ll in seq(along=colorMaps)) {
      colorMap <- colorMaps[[ll]];
      tags <- colorMap$tags;
      verbose && enter(verbose, sprintf("Color map #%d ('%s')", ll, tags));
#      verbose && str(verbose, colorMap$transform);
#      verbose && str(verbose, colorMap$palette);
      writeImage(df, path=path, transform=colorMap$transform, 
                                   palette=colorMap$palette, tags=tags, ...);
      verbose && exit(verbose);
    }
    verbose && exit(verbose);
  }

  invisible(path);
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
setMethodS3("process", "ArrayExplorer", function(this, arrays=NULL, ..., verbose=FALSE) {
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

  
  verbose && enter(verbose, "Generating ", class(this)[1], " report");

  # Setup HTML, CSS, Javascript files first
  setup(this, ..., verbose=less(verbose));

  # Generate bitmap images
  writeImages(this, arrays=arrays, ..., verbose=less(verbose));

  # Update samples.js
  updateSamplesFile(this, ..., verbose=less(verbose));

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
setMethodS3("display", "ArrayExplorer", function(this, ..., verbose=FALSE) {
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
  path <- getPath(this);
  path <- getParent(path);
  path <- getParent(path);
  filename <- sprintf("%s.html", class(this)[1]);
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
# 2007-02-08
# o Created from ChromosomeExplorer.R.
##############################################################################
