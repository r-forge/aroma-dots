###########################################################################/**
# @RdocClass SpatialReporter
#
# @title "The SpatialReporter class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixCelSetReporter".}
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
setConstructorS3("SpatialReporter", function(...) {
  extend(AffymetrixCelSetReporter(...), "SpatialReporter"
  )
})

setMethodS3("as.character", "SpatialReporter", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, paste("Number of arrays:", nbrOfArrays(this)));
  colorMaps <- getColorMaps(this);
  if (length(colorMaps) == 0)
    colorMaps <- "<no color maps; set before processing>";
  s <- c(s, paste("Color maps:", paste(colorMaps, collapse="; ")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getReportSet", "SpatialReporter", function(this, ...) {
  "spatial";
}, protected=TRUE)


setMethodS3("addColorMap", "SpatialReporter", function(this, colorMap, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'colorMap':
  colorMap <- Arguments$getCharacter(colorMap, nchar=c(1,Inf));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Parser argument 'colorMap'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  parts <- strsplit(colorMap, split=",")[[1]];
  tags <- paste(parts, collapse=",");
  n <- length(parts);
  if (n > 3) {
    throw("Argument 'colorMap' must not contain more than three parts: ", 
                                                                   tags);
  }
  if (n == 1) {
    transforms <- list();
  } else {
    transforms <- as.list(parts[-n]);
  }
  palette <- parts[n];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Setup transforms
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (length(transforms) == 0) {
    transforms <- list(sqrt);
  } else {
    transforms <- lapply(transforms, FUN=function(transform) {
      # Check transform
      if (!exists(transform, mode="function")) {
        throw("Argument 'colorMap' specifies an unknown transform function ('",
                                                 transform, "'): ", colorMap);
      }
      get(transform, mode="function");
    })
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Setup palette
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (is.null(palette)) {
    palette <- gray.colors(256);
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
    transforms = transforms,
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


setMethodS3("setColorMaps", "SpatialReporter", function(this, colorMaps=c("sqrt,yellow", "sqrt,rainbow"), ...) {
  this$.colorMaps <- NULL;
  for (colorMap in colorMaps) {
    addColorMap(this, colorMap, ...);
  }
})

setMethodS3("getColorMaps", "SpatialReporter", function(this, parsed=FALSE, ...) {
  colorMaps <- this$.colorMaps;
  if (!parsed)
    colorMaps <- sapply(colorMaps, .subset2, "tags");
  colorMaps;
})


setMethodS3("writeImages", "SpatialReporter", function(this, aliases=NULL, ...) {
  # Get the CEL set of interest
  cs <- getDataSet(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(aliases)) {
    aliases <- Arguments$getCharacters(aliases, length=nbrOfArrays(cs));
  }
  
  # Get the path to the image directory
  path <- getPath(this);

  # Get the color maps to be generated
  colorMaps <- getColorMaps(this, parsed=TRUE);
  if (length(colorMaps) == 0) {
    warning("No color maps specified. Nothing to do.");
    return(invisible(path));
  }

  # For each array...
  for (kk in seq(cs)) {
    df <- getFile(cs, kk);
    setAlias(df, aliases[kk]);
    verbose && enter(verbose, sprintf("Array #%d ('%s')", kk, getName(df)));
    verbose && cat(verbose, "Alias: ", getAlias(df));
    # For each color map...
    for (ll in seq(along=colorMaps)) {
      colorMap <- colorMaps[[ll]];
      tags <- colorMap$tags;
      verbose && enter(verbose, sprintf("Color map #%d ('%s')", ll, tags));
#      verbose && str(verbose, colorMap$transforms);
#      verbose && str(verbose, colorMap$palette);
      writeImage(df, path=path, transforms=colorMap$transforms, 
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
setMethodS3("process", "SpatialReporter", function(this, ..., verbose=FALSE) {
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

  # Generate bitmap images
  writeImages(this, ..., verbose=less(verbose));

  verbose && exit(verbose);
})



##############################################################################
# HISTORY:
# 2007-03-19
# o Updated addColorMap() to accept multiple transforms.
# o Created from ArrayExplorer.R.
##############################################################################
