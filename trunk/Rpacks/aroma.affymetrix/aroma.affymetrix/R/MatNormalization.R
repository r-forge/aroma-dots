###########################################################################/**
# @RdocClass MatNormalization
#
# @title "The MatNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents a normalization method that corrects for systematic
#  effects in the probe intensities due to differences in the number of
#  A, C, G, and T:s and the match scores according to MAT - model-based
#  analysis of tiling arrays Johnson et al. PNAS 2006
#  (senior author: Shirley Liu, Harvard).
# }
# 
# @synopsis 
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of 
#     @see "AbstractProbeSequenceNormalization".}
#   \item{model}{A @character string specifying the model used to fit 
#     the base-count effects.}
#   \item{numChunks}{The number of chunks to split the data into to fit the model}
#   \item{numBins}{The number of bins to use for the variance smoothing step}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \section{Requirements}{
#   This class requires that an aroma probe sequence file and aroma
#   match scores file is available for the chip type.
# }
# 
# @author
#*/###########################################################################
setConstructorS3("MatNormalization", function(..., unitsToFit=NULL, model=c("lm"), numChunks=25, numBins=200) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'model':
  model <- match.arg(model);

  extend(AbstractProbeSequenceNormalization(..., unitsToFit=unitsToFit), "MatNormalization",
    .model = model,
    .scaleResiduals = TRUE,
    .numChunks = as.integer(numChunks),
    .numBins = as.integer(numBins)
  )
})


setMethodS3("getAromaCellMatchscoreFile", "MatNormalization", function(this, ..., force=FALSE) {
  apm <- this$.apm;

  if (force || is.null(apm)) {
    dataSet <- getInputDataSet(this);
    cdf <- getCdf(dataSet);
    chipType <- getChipType(cdf, fullname=FALSE);
    apm <- AromaCellMatchscoreFile$byChipType(chipType, ...);
    this$.apm <- apm;
  }

  apm;
}, protected=TRUE)



setMethodS3("getAsteriskTags", "MatNormalization", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", this, collapse=collapse, ...);

  # Add model tag?
  model <- this$.model;
  tags <- c(tags, model);

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, private=TRUE)


setMethodS3("getParameters", "MatNormalization", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  params <- c(params, list(
    model = this$.model,
    numChunks = this$.numChunks
  ));

  params;
}, private=TRUE)



setMethodS3("getDesignMatrix", "MatNormalization", function(this, cells=NULL, model=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  if (is.null(cells)) {
  } else {
    # Validated below...
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving design matrix");
  verbose && cat(verbose, "Cells:");
  verbose && str(verbose, cells);

  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Locating match scores annotation data");
  # Locate AromaCellMatchscoreFile holding match scores
  apm <- getAromaCellMatchscoreFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Read sequence matrix");
  sm<-readSequenceMatrix(aps, cells=cells, verbose=verbose);
  verbose && exit(verbose);
  verbose && enter(verbose, "Read match scores");
  ms<-readColumns(apm, rows=cells, verbose=verbose);
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Construct design matrix");
  nT <- rowSums(sm == "T");
  G <- (sm == "G")+0;
  A <- (sm == "A")+0;
  C <- (sm == "C")+0;
  designMatrix <- cbind(nT, A, C, G, rowSums(A)^2, rowSums(C)^2, rowSums(G)^2, nT^2, log(as.integer(ms[,1])));
  
  # Garbage collect
  rm(nT,G,A,C,ms,sm);
  gc <- gc();
  verbose && print(verbose, gc);
  
  verbose && str(verbose, designMatrix);
  verbose && cat(verbose, "object.size(designMatrix): ", object.size(designMatrix));
  verbose && exit(verbose);

  verbose && exit(verbose);

  designMatrix;
}, private=TRUE)



setMethodS3("fitOne", "MatNormalization", function(this, df, ..., verbose=FALSE) {

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Locating match scores annotation data");
  # Locate AromaCellMatchscoreFile holding match scores
  apm <- getAromaCellMatchscoreFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading 'non-missing' cells to fit");
  cellsToFit<-whichVector( !(isMissing(aps) | isMissing(apm)) );
  verbose && cat(verbose, "Cells to fit:");
  verbose && str(verbose, cellsToFit);
  verbose && exit(verbose);

  numChunks<-this$.numChunks

  # this code adopted from Dave Fourniers 17/08/2007 post to r-help mailing list
  # entitled "[R] Linear models over large datasets"

  incr<-ceiling(length(cellsToFit)/numChunks)+1
  
  verbose && enter(verbose, "Reading signals to fit");
  y <- extractMatrix(df, cells=cellsToFit, verbose=less(verbose, 10));
  verbose && exit(verbose);

  verbose && enter(verbose, "Log2 transforming signals");
  y <- log2(y);
  verbose && cat(verbose, "Target log2 probe signals:");
  verbose && str(verbose, y);
  verbose && exit(verbose);
 
  start <- xtx <- xty <- 0
   
  while(start < length(cellsToFit)) {
  
    verbose && enter(verbose, "Working on indices over range");
    indSubset <- seq(start + 1, min(start + incr, length(cellsToFit)));
    rng <- range(indSubset)
    verbose && cat(verbose, paste(rng/length(cellsToFit),sep=""));
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading design matrix for this iteration");
    X<-getDesignMatrix(this, cells=cellsToFit[indSubset], verbose=verbose);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating cross products");
    xtx <- xtx + crossprod(X);
    xty <- xty + crossprod(X, y[indSubset]);
    verbose && exit(verbose);

    start <- start + incr;
  }
  
  verbose && enter(verbose, "Solving linear system");
  #fit <- list(xtx=xtx, xty=xty) #,beta=solve(xtx, xty))
  verbose && exit(verbose);
  fit <- list(beta=solve(xtx, xty), scaleResiduals=this$.scaleResiduals);

  fit;
}, protected=TRUE)


setMethodS3("predictOne", "MatNormalization", function(this, fit, ..., verbose=FALSE) {

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Locating match scores annotation data");
  # Locate AromaCellMatchscoreFile holding match scores
  apm <- getAromaCellMatchscoreFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Allocating mu vector");
  mu <- double(nbrOfCells(aps));
  verbose && str(verbose, mu);
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Reading 'non-missing' cells to predict");
  cellsToPredict <- whichVector( !(isMissing(aps) | isMissing(apm)) );
  verbose && cat(verbose, "Cells to predict:");
  verbose && str(verbose, cellsToPredict);
  verbose && exit(verbose);

  numChunks<-this$.numChunks

  incr<-ceiling(length(cellsToPredict)/numChunks)+1
  start <- 0
     
  while(start < length(cellsToPredict)) {
  
    verbose && enter(verbose, "Working on indices over range");
    indSubset <- seq(start + 1, min(start + incr, length(cellsToPredict)));
    rng <- range(indSubset);
    # HB: cat(str(paste(...))) ?!?
    verbose && cat(verbose, str(paste(round(rng/length(cellsToPredict),3), sep="")));
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading design matrix for this iteration");
    X <- getDesignMatrix(this, cells=cellsToPredict[indSubset], verbose=verbose)
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating mu");
    mu[ cellsToPredict[indSubset] ] <- X %*% fit$beta;
    verbose && exit(verbose);

    start <- start + incr;
  }
  
  rm(X,indSubset);
  gc <- gc();
  
  #y <-
  
  #numBins <- this$.numBins
  #q <- quantile(mu[cellsToPredict],prob=(0:numBins)/numBins)
  #cuts<-cut(mu[cellsToPredict],breaks=q,labels=1:(length(q)-1))  # define
  #ss<-split(data.frame(resid),cuts)
  #ssvar<-sapply(ss,var)
  #v<-ssvar[as.character(cuts)]
  #for(j in 1:length(b))
  #rr<-resid/sqrt(v)

  mu;
}, protected=TRUE)

###########################################################################/**
# @RdocMethod process
#
# @title "Normalizes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already normalized is re-normalized, 
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "MatNormalization", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Normalization data set for probe-sequence effects");

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
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get (and create) the output path
  outputPath <- getPath(this);

  verbose && enter(verbose, "Locating probe-sequence annotation data");
  # Locate AromaCellSequenceFile holding probe sequences
  aps <- getAromaCellSequenceFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Locating match scores annotation data");
  # Locate AromaCellMatchscoreFile holding match scores
  apm <- getAromaCellMatchscoreFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading 'non-missing' cells to fit");
  cellsToFit <- whichVector( !(isMissing(aps) | isMissing(apm)) );
  verbose && cat(verbose, "Cells to fit:");
  verbose && str(verbose, cellsToFit);
  verbose && exit(verbose);

  numChunks <- this$.numChunks;
  numBins <- this$.numBins;
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize all arrays simultaneously
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfArrays <- nbrOfArrays(ds);
  df <- getFile(ds, 1);
  nbrOfCells <- nbrOfCells(df);
  verbose && enter(verbose, "Normalizing ", nbrOfArrays, " arrays");
  verbose && enter(verbose, "Path: ", outputPath);
  
  start <- xtx <- 0;
  xty <- lapply(seq_len(nbrOfArrays), FUN=function(u) return(0));
  incr <- ceiling(length(cellsToFit)/numChunks)+1;
   
  while(start < length(cellsToFit)) {
    verbose && enter(verbose, "Working on indices over range");
    indSubset <- seq(start + 1, min(start + incr, length(cellsToFit)))
    rng <- range(indSubset);
    verbose && cat(verbose, paste(rng/length(cellsToFit), sep=""));
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading design matrix and data for this iteration");
    X <- getDesignMatrix(this, cells=cellsToFit[indSubset], verbose=verbose);
    Y <- extractMatrix(ds, cells=cellsToFit[indSubset], verbose=verbose);
    verbose && exit(verbose);

    verbose && enter(verbose, "Log2 transforming signals");
    Y <- log2(Y);
    verbose && cat(verbose, "Target log2 probe signals:");
    verbose && str(verbose, Y);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating cross products");
    xtx <- xtx + crossprod(X);
    for (kk in seq_len(nbrOfArrays)) {
      xty[[kk]] <- xty[[kk]] + crossprod(X, Y[,kk]);
    }
    verbose && exit(verbose);

    start <- start + incr;
  }  # while (start < ...)

  fits <- lapply(xty, FUN=function(u) list(beta=solve(xtx,u)));
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save model fits
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in seq_len(nbrOfArrays)) {
    df <- getFile(ds, kk);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                              kk, getName(df), nbrOfArrays));

    fullname <- getFullName(df);

    # Store model fit 
    verbose && enter(verbose, "Saving model fit");
    # Store fit and parameters (in case someone are interested in looking
    # at them later; no promises of backward compatibility though).
    filename <- sprintf("%s,fit.RData", fullname);
    fitPathname <- Arguments$getWritablePathname(filename, 
                                                    path=outputPath, ...);
    saveObject(fits[[kk]], file=fitPathname);
    verbose && str(verbose, fits[[kk]], level=-50);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);
  } # for (kk ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate residuals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  start <- xtx <- 0;
  mu <- vector("list", nbrOfArrays);
  incr <- ceiling(length(cellsToFit)/numChunks)+1;
   
  while(start < length(cellsToFit)) {
    verbose && enter(verbose, "Working on indices over range");
    indSubset <- seq(start + 1, min(start + incr, length(cellsToFit)));
    rng <- range(indSubset);
    verbose && cat(verbose, paste(rng/length(cellsToFit),sep=""));
    verbose && exit(verbose);

    verbose && enter(verbose, "Reading design matrix for this iteration");
    X <- getDesignMatrix(this, cells=cellsToFit[indSubset], verbose=verbose);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating model fits");
    xtx <- xtx + crossprod(X);
    for (kk in seq_len(nbrOfArrays)) {
       mu <- X %*% fits[[kk]]$beta;
       verbose && str(verbose, mu);

       df <- getFile(ds, kk);
       verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                              kk, getName(df), nbrOfArrays));

      fullname <- getFullName(df);
      filename <- sprintf("%s.CEL", fullname);
      pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");
      createFrom(df, filename=pathname, path=NULL, verbose=less(verbose));
      verbose && exit(verbose);

      updateCel(pathname, indices=cellsToFit[indSubset], intensities=2^as.numeric(mu), verbose=TRUE);

      verbose && exit(verbose);
    } # for (kk ...)

    rm(mu);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    verbose && exit(verbose);

    start <- start + incr;
    #start <- length(cellsToFit)+1;

  } # while (start < ...)
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Scale residuals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in seq_len(nbrOfArrays)) {
  
      verbose && enter(verbose, "Binning predicted values, calculating and scaling residuals");
      df <- getFile(ds, kk);

      y <- log2(extractMatrix(df, cells=cellsToFit, verbose=verbose))

      fullname <- getFullName(df);
      filename <- sprintf("%s.CEL", fullname);
      pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

      mu <- log2(readCel(pathname, indices=cellsToFit, readOutliers=FALSE, readHeader=FALSE, readMasked=FALSE, verbose=less(verbose,10))$intensities);
      r <- y - mu;

      q <- quantile(mu, prob=(0:numBins)/numBins);
      cuts <- cut(mu,breaks=q,labels=1:(length(q)-1));  # define
      ss <- split(r, cuts);
      ssvar <- sapply(ss, var);
      v <- ssvar[as.character(cuts)];
      r <- r/sqrt(v);

      #return(list(y=y,mu=mu,r=r))

      # Create CEL file to store results, if missing
      verbose && enter(verbose, "Creating CEL file for results, if missing");
      createFrom(df, filename=pathname, path=NULL, verbose=less(verbose));
      verbose && exit(verbose);

      updateCel(pathname, indices=cellsToFit, intensities=2^as.numeric(r), verbose=TRUE);
      rm(q,ss,ssvar,v,r,y);
      gc <- gc();
      verbose && print(verbose, gc);

      verbose && exit(verbose);
  } # for (kk ...)

  outputDataSet <- getOutputDataSet(this, force=TRUE);

  verbose && exit(verbose);
  
  invisible(outputDataSet);
})






############################################################################
# HISTORY:
# 2008-10-29 [MR]
# o Created from BaseCountNormalization.R
############################################################################
