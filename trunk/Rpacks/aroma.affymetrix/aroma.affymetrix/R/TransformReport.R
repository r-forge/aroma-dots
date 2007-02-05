###########################################################################/**
# @RdocClass TransformReport
#
# @title "The TransformReport class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis 
#
# \arguments{
#   \item{inSet}{The input data set as an @see "AffymetrixCelSet".}
#   \item{outSet}{The output data set as an @see "AffymetrixCelSet".}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \details{
# }
#
# @author
#*/###########################################################################
setConstructorS3("TransformReport", function(inSet=NULL, outSet=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'inSet':
  if (!is.null(inSet)) {
    if (!inherits(inSet, "AffymetrixCelSet")) {
      throw("Argument 'inSet' is not an AffymetrixCelSet object: ", 
                                                            class(inSet)[1]);
    }

    if (!inherits(outSet, "AffymetrixCelSet")) {
      throw("Argument 'outSet' is not an AffymetrixCelSet object: ", 
                                                           class(outSet)[1]);
    }

    # Check for compatibility
#    if (!equals(getCdf(inSet), getCdf(outSet))) {
#      throw("Argument 'inSet' and 'outSet' have incompatible CDFs.");
#    }
  }

  extend(Object(), "TransformReport", 
    .inSet = inSet,
    .outSet = outSet
  )
}, abstract=TRUE)


setMethodS3("clearCache", "TransformReport", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c()) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
}, private=TRUE)


setMethodS3("getRootPath", "TransformReport", function(this, ...) {
  sprintf("pp%s", capitalize(class(this)[1]));
}, private=TRUE)


setMethodS3("as.character", "TransformReport", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  ds <- getInputDataSet(this);
  s <- c(s, sprintf("Input data set: %s", getFullName(ds)));
  ds <- getOutputDataSet(this);
  s <- c(s, sprintf("Output data set: %s", getFullName(ds)));
  s <- c(s, sprintf("Number of arrays: %d (%.2fMB)", 
                           nbrOfArrays(ds), getFileSize(ds)/1024^2));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(ds))));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the output data set"
#
# \description{
#  @get "title", which is the same as the input data set.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "TransformReport", function(this, ...) {
  ds <- getOutputDataSet(this);
  getName(ds);
})


###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the output data set"
#
# \description{
#  @get "title", which equals the tags of the input data set plus the tags
#  of this transformation.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
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
setMethodS3("getTags", "TransformReport", function(this, ...) {
  tags <- this$.tags;

  ds <- getOutputDataSet(this);
  tags <- getTags(ds);

  tags;
})



###########################################################################/**
# @RdocMethod getFullName
#
# @title "Gets the full name of the output data set"
#
# \description{
#  @get "title", which is the name with comma separated tags.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFullName", "TransformReport", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path of the output data set"
#
# \description{
#  @get "title".
#  If non-existing, then the directory is created.
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
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPath", "TransformReport", function(this, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  ds <- getOutputDataSet(this);
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf);
  chipType <- gsub("[,-]monocell$", "", chipType);  # AD HOC? /HB 2006-12-08

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
})


###########################################################################/**
# @RdocMethod getInputDataSet
#
# @title "Gets the source data set"
#
# \description{
#  @get "title" that is to be (or has been) transformed.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getInputDataSet", "TransformReport", function(this, ...) {
  this$.inSet;
})



###########################################################################/**
# @RdocMethod getOutputDataSet
#
# @title "Gets the transformed data set"
#
# \description{
#  @get "title", if processed.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCelSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getOutputDataSet", "TransformReport", function(this, ...) { 
  this$.outSet;
})


setMethodS3("getYY", "TransformReport", function(this, array, transform=NULL, subset=1/8, field="intensities", ...) {
  inSet <- getInputDataSet(this);

  cdf <- getCdf(inSet);
  if (length(subset) == 1) {
    indices <- seq(from=1, to=nbrOfCells(cdf), length=subset*nbrOfCells(cdf));
    indices <- as.integer(indices);
  } else if (length(subset) > 1) {
    indices <- Arguments$getIndices(subset, range=c(1, nbrOfCells(cdf)));
  } else {
    indices <- subset;
  }

  outSet <- getOutputDataSet(this);
  df1 <- getFile(inSet, array);
  df2 <- getFile(outSet, array);
  res <- list(
    y1 = getData(df1, indices=indices, ..., fields=field)[[field]],
    y2 = getData(df2, indices=indices, ..., fields=field)[[field]],
    array = array,
    df1 = df1,
    df2 = df2
  );

  if (!is.null(transform)) {
    res$y1 <- transform(res$y1);
    res$y2 <- transform(res$y2);
  }

  res;
})



setMethodS3("plotYYSpline", "TransformReport", function(this, xlim=c(0,65535), xlab=expression(y), ylab=expression(tilde(y)==h(y)), dcol="#cccccc", main=NULL, ...) {
  suppressWarnings({
    yy <- getYY(this, ...);
  })

  if (is.null(main)) {
    main <- getName(yy$df2);
  }


  suppressWarnings({
    fit <- plotXYSpline(yy$y1, yy$y2, lwd=4, dcol=dcol, xlim=xlim, xlab=xlab, ylab=ylab, ...);
  })

  cdf <- getCdf(yy$df2);
  stextChipType(cdf, line=-1);
  stextSize(yy$df2, size=length(yy$y1));

  invisible(fit);
})

setMethodS3("plotYYSplineLog2", "TransformReport", function(this, xlim=c(0,16), xlab=expression(log[2](y)), ylab=expression(log[2]*tilde(y)==log[2]*h(y)), ...) {
  plotYYSpline(this, transform=log2, xlim=xlim, xlab=xlab, ylab=ylab, ...);
})



setMethodS3("writeImages", "TransformReport", function(this, path=NULL, width=640, height=width, ..., skip=TRUE, verbose=FALSE) {
  pngDev <- System$findGraphicsDevice();
 
  rootPath <- getRootPath(this);
  name <- getName(this);
  tags <- getTags(this);
  tags <- paste(tags, collapse=",");
  path <- file.path(rootPath, name, tags);
  path <- Arguments$getWritablePath(path);

  outSet <- getOutputDataSet(this);
  nbrOfArrays <- nbrOfArrays(outSet);

  verbose && enter(verbose, "Writing images for ", nbrOfArrays, " arrays");

  verbose && printf(verbose, "Image dimension: %.0fx%.0f\n", width, height);

  for (kk in seq(length=nbrOfArrays)) {
    df <- getFile(outSet, kk);
    fullname <- getFullName(df);
    verbose && enter(verbose, "Output CEL file: ", fullname);

    # Plot (y,y)
    tags <- c("YvY");
    imgname <- paste(c(fullname, tags), collapse=",");
    filename <- sprintf("%s.png", imgname);
    pathname <- file.path(path, filename);

    verbose && cat(verbose, "Image pathname: ", pathname);
    if (!skip || !isFile(pathname)) {
      pngDev(pathname, width=width, height=height);
      tryCatch({
        plotYYSpline(this, array=kk, ...);
      }, finally = {
        dev.off();
      })
    }

    # Plot (log2(y),log2(y))
    tags <- c("YvY,log2");
    imgname <- paste(c(fullname, tags), collapse=",");
    filename <- sprintf("%s.png", imgname);
    pathname <- file.path(path, filename);

    verbose && cat(verbose, "Image pathname: ", pathname);
    if (!skip || !isFile(pathname)) {
      pngDev(pathname, width=width, height=height);
      tryCatch({
        plotYYSplineLog2(this, array=kk, ...);
      }, finally = {
        dev.off();
      })
    }

    # Garbage collection
    if (kk %% 10 == 0)
      verbose && print(verbose, gc());

    verbose && exit(verbose);
  }

  # Garbage collection
  verbose && print(verbose, gc());

  verbose && exit(verbose);
})




############################################################################
# HISTORY:
# 2007-02-04
# o Created.
############################################################################
