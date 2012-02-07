setConstructorS3("IlluminaIDatSet", function(...) {
  extend(GenericDataFileSet(...), "IlluminaIDatSet");
})

setMethodS3("getPlatform", "IlluminaIDatSet", function(static, ...) {
  IlluminaIDatFile$getPlatform(...);
}, static=TRUE);


setMethodS3("extractMatrix", "IlluminaIDatSet", function(this, drop=FALSE, ...) {
  # Argument 'drop':
  drop <- Arguments$getLogical(drop);

  sampleNames <- getNames(this);
  nbrOfSamples <- length(sampleNames);

  Y <- NULL;
  for (ii in seq(length=nbrOfSamples)) {
    df <- getFile(this, ii);
    y <- extractMatrix(df, ...);
    unitNames <- rownames(y);
    if (ii == 1L) {
      dimnames <- list(unitNames, sampleNames);
      naValue <- NA;
      storage.mode(naValue) <- storage.mode(y);
      Y <- matrix(naValue, nrow=length(y), ncol=nbrOfSamples, 
                           byrow=FALSE, dimnames=dimnames);
    } else {
      # Sanity check
      stopifnot(identical(unitNames, rownames(Y)));
    }
    Y[,ii] <- y[,1L];
    rm(y);
  } # for (ii ...)

  if (drop && (nbrOfSamples == 1)) {
    Y <- Y[,1L, drop=TRUE];
  }

  Y;
}, protected=TRUE)



setConstructorS3("IlluminaIDatFile", function(...) {
  extend(GenericDataFile(...), "IlluminaIDatFile");
})

setMethodS3("getPlatform", "IlluminaIDatFile", function(static, ...) {
  "Illumina";
}, static=TRUE);


setMethodS3("readQuants", "IlluminaIDatFile", function(this, ...) {
  pathname <- getPathname(this);
  data <- readIDAT(pathname, ...);
  data <- data$Quants;
  data;
}, protected=TRUE)


setMethodS3("readSignals", "IlluminaIDatFile", function(this, ...) {
  data <- readQuants(this, ...);
  data <- data[,"Mean"];
  data;
}, protected=TRUE)


setMethodS3("extractMatrix", "IlluminaIDatFile", function(this, drop=FALSE, ...) {
  # Argument 'drop':
  drop <- Arguments$getLogical(drop);

  y <- readSignals(this, ...);

  if (!drop) {
    y <- as.matrix(y);
  }
  y;
}, protected=TRUE)




############################################################################
# HISTORY:
# 2012-02-07
# o Created.
############################################################################
