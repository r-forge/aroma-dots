setConstructorS3("SpatialRowColumnNormalization", function(..., spar=c(0.7,0.7), h=c(20,20)) {
  # Argument 'spar':
  spar <- Arguments$getDoubles(spar, range=c(0,Inf));

  # Argument 'h':
  h <- Arguments$getIntegers(h, range=c(1,Inf));


  extend(ProbeLevelTransform(...), "SpatialRowColumnNormalization",
    .refFile = NULL,
    .refData = NULL,
    .spar = spar,
    .h = h
  );
})


setMethodS3("getParameters", "SpatialRowColumnNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  # Get parameters of this class
  params2 <- list(
    spar = getSpar(this),
    h = getH(this),
    yR = getReferenceData(this)
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, private=TRUE)
 

setMethodS3("getSpar", "SpatialRowColumnNormalization", function(this, ...) {
  this$.spar;
})

setMethodS3("getH", "SpatialRowColumnNormalization", function(this, ...) {
  this$.h;
})

setMethodS3("getReferenceFile", "SpatialRowColumnNormalization", function(this, force=FALSE, ...) {
  refFile <- this$.refFile;

  if (force || is.null(refFile)) {
    ds <- getInputDataSet(this);
    refFile <- getAverageFile(ds, ...);
    this$.refFile <- refFile;
  }

  refFile;
})


setMethodS3("getReferenceData", "SpatialRowColumnNormalization", function(this, force=FALSE, ...) {
  refData <- this$.refData;
  if (force || is.null(refData)) {
    refFile <- getReferenceFile(this, ...);
    refData <- readRawDataRectangle(refFile, field="intensities", drop=TRUE, ...);
    this$.refData <- refData;
  }
  refData;
})



setMethodS3("process", "SpatialRowColumnNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Normalizing data set spatially in blocks of rows and columns");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already normalized");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve/calculate the reference file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving reference file");
  refFile <- getReferenceFile(this, verbose=less(verbose));
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);

  # Get algorithm parameters (including the target distribution above)
  params <- getParameters(this);
  yR <- params$yR;

  verbose && cat(verbose, "Cell indices:");
  cells <- matrix(seq(along=yR), nrow=nrow(yR), ncol=ncol(yR), byrow=TRUE);
  verbose && str(verbose, cells);

  # Get the output path
  outputPath <- getPath(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- nbrOfArrays(ds);
  verbose && enter(verbose, "Normalizing ", nbrOfArrays, " arrays");
  verbose && enter(verbose, "Path: ", outputPath);
  dataFiles <- list();
  for (kk in seq_len(nbrOfArrays)) {
    df <- getFile(ds, kk);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                              kk, getName(df), nbrOfArrays));

    fullname <- getFullName(df);
    filename <- sprintf("%s.CEL", fullname);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...); 

    # Already normalized?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Normalized data file already exists: ", pathname);
    } else {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Reading data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Reading all probe intensities");
      y <- readRawDataRectangle(df, field="intensities", drop=TRUE, ...);
      verbose && str(verbose, y);
      verbose && exit(verbose);

      verbose && enter(verbose, "Transforming to the log-scale");
      y <- log2(y);
      verbose && str(verbose, y);
      verbose && exit(verbose);


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Normalizing log-ratios data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Normalizing rows and columns in blocks");
      y <- norm2d(y, yTarget=yR, spar=params$spar, h=params$h, ...);
      verbose && str(verbose, y); 
      verbose && exit(verbose); 

      verbose && enter(verbose, "Back-transforming to intensity scale");
      y <- 2^y;
      y <- as.vector(y);
      verbose && str(verbose, y); 
      verbose && exit(verbose); 

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Storing data
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Storing normalized data");
    
      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");
      createFrom(df, filename=pathname, path=NULL, verbose=less(verbose));
      verbose && exit(verbose);

      # Write calibrated data to file
      verbose2 <- -as.integer(verbose)-2;
      verbose && str(verbose, cells);
      updateCel(pathname, indices=cells, intensities=y, verbose=verbose2);
      rm(y, verbose2);

      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose); 
    } # if-else

    # Retrieving normalized data file
    dfN <- newInstance(df, pathname);

    # CDF inheritance
    setCdf(dfN, cdf);

    # Record 
    dataFiles[[kk]] <- dfN;

    rm(df, dfN);
    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  # Garbage collect
  rm(dataFiles, ds);
  gc <- gc();
  verbose && print(verbose, gc); 

  # Update the output data set
  outputDataSet <- getOutputDataSet(this, force=TRUE);

  verbose && exit(verbose);
  
  invisible(outputDataSet);
}) 


############################################################################
# HISTORY: 
# 2008-04-02
# o Note: Normalizing log-ratios and then transforming back to chip effects
#   might not do.  See my PPT notes from today.
# 2008-03-19
# o Created.
############################################################################ 
