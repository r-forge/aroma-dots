setConstructorS3("GenePixResultSet", function(...) {
  extend(GenericDataFileSet(...), "GenePixResultSet");
})



setMethodS3("findByName", "GenePixResultSet", function(static, ..., chipType=NULL, paths=c("rawData(|,.*)/")) {
  # Arguments 'chipType':`
  if (!is.null(chipType)) {
    chipType <- Arguments$getCharacter(chipType);
  }

  # Arguments 'paths':
  if (is.null(paths)) {
    paths <- eval(formals(findByName.GenePixResultSet)[["paths"]]);
  }

  # NextMethod() does not work here.
  findByName.GenericDataFileSet(static, ..., subdirs=chipType, paths=paths);
}, static=TRUE)



setMethodS3("byName", "GenePixResultSet", function(static, name, tags=NULL, chipType=NULL, paths=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'chipType':
  if (!is.null(chipType)) {
    chipType <- Arguments$getCharacter(chipType);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Setting up ", class(static)[1], " by name");

  verbose && cat(verbose, "Name: ", name);
  verbose && cat(verbose, "Tags: ", paste(tags, collapse=","));

  suppressWarnings({
    paths <- findByName(static, name, tags=tags, chipType=chipType, 
                       paths=paths, firstOnly=FALSE, ...);
  })
  if (is.null(paths)) {
    path <- file.path(paste(c(name, tags), collapse=","), chipType);
    throw("Cannot create ", class(static)[1], ".  No such directory: ", path);
  }

  verbose && cat(verbose, "Paths to possible data sets:");
  verbose && print(verbose, paths);

  res <- NULL;
  for (kk in seq(along=paths)) {
    path <- paths[kk];
    verbose && enter(verbose, sprintf("Trying path #%d of %d", kk, length(paths)));
    verbose && cat(verbose, "Path: ", path);

    tryCatch({
      suppressWarnings({
        res <- byPath(static, path=path, ..., verbose=verbose);
      });
    }, error = function(ex) {
      verbose && cat(verbose, "Data set could not be setup for this path, because:");
      verbose && cat(verbose, ex$message);
    });

    if (!is.null(res)) {
      if (nbrOfFiles(res) > 0) {
        verbose && cat(verbose, "Successful setup of data set.");
        verbose && exit(verbose);
        break;
      }
    }

    verbose && exit(verbose);
  } # for (kk ...)

  if (is.null(res)) {
    throw(sprintf("Failed to setup a data set for any of %d data directories located.", length(paths)));
  }

  verbose && exit(verbose);

  res;
}, static=TRUE) # byName()




setMethodS3("byPath", "GenePixResultSet", function(static, ..., pattern=NULL, fileClass=getFileClass(static)) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fileClass':
  clazz <- Class$forName(fileClass);
  dfStatic <- getStaticInstance(clazz);
  dfStatic <- Arguments$getInstanceOf(dfStatic, getFileClass(static), .name="fileClass");

  # Argument 'pattern':
  if (!is.null(pattern)) {
    pattern <- Arguments$getRegularExpression(pattern);
  }

  # Default filename extension pattern
  if (is.null(pattern)) {
    pattern <- getExtensionPattern(dfStatic);
  }

  byPath.GenericDataFileSet(static, ..., pattern=pattern, fileClass=fileClass);
}, static=TRUE);


setMethodS3("getPlatform", "GenePixResultSet", function(this, ...) {
  df <- getFile(this, 1);
  getPlatform(df);
})

setMethodS3("getChipType", "GenePixResultSet", function(this, ...) {
  df <- getFile(this, 1);
  getChipType(df);
})


############################################################################
# HISTORY:
# 2011-05-23
# o Added findByName() and byName().
# o Added readLimmaData() from old GenePixDataSet.R.
# o Created.
############################################################################
