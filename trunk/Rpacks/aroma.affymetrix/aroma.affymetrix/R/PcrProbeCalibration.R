###########################################################################/**
# @RdocClass PcrProbeCalibration
#
# @title "The PcrProbeCalibration class"
#
# \description{
#  @classhierarchy
#
#  This class represents the first step of pre-processing that corrects for
#  PCR effects in the \emph{probe-level} intensities due to probe sequences 
#  and PCR fragment lengths, cf. Carvalho et al (2006).
# }
# 
# @synopsis 
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of
#      @see "ProbePreprocessing".}
#   \item{subsetToFit}{The fraction of PM probes to be used for fitting
#     the calibration function.}
#   \item{flavor}{A @character string specifying what method to immitate.
#     If \code{"oligo"}, the calibration gives identical results to the
#     \pkg{oligo} package; the \code{subsetToFit} argument is ignored.
#   }
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
#
# \details{
#   For sample \eqn{i=1,...,I} and probe \eqn{k=1,...,K}, the PM intensity
#   is modelled as:
#   \deqn{
#     \log_2(PM_{ik}) = \mu_{ik} + g(L_k) + \sum_{t=1}^{25} 
#                       \sum_{b} h_b(t) I(b_{kt}=b) + \xi_{ik},
#   }
#   where  \eqn{\xi} is \eqn{N(0,\sigma^2)}. 
#   The \eqn{g(L_k)} is the PCR fragment-length effect for a probe in a 
#   SNP of length \eqn{L_k}.  Thus, all probes in a SNP have the same
#   fragment-length effect \eqn{g(L_k)} (independ on sample).
#   The \eqn{h_b(t)} is the effect of having base \eqn{b=\{A,C,G,T\}}
#   at nucleotide position \eqn{t=1,...,25}.
#   The \eqn{\mu_{ik}} is the signal of interest.
# }
#
# \section{Relationship to the oligo package}{
#  The methods of class are adopted from the @see "oligo::snprma" function 
#  in the \pkg{oligo} package. By setting argument 
#  \code{subsetToFit="asOligo"}, this calibration method produces the same
#  result as the \pkg{oligo} package.
# }
#
# \section{Requirments}{
#   To use this class the platform-design package for the chip type of the
#   data set as well as the \pkg{oligo} package must be is installed.
#   The methods, which are adopted from the \pkg{oligo} package, have been
#   rewritten to minimize the memory usage. For instance, the platform-design
#   package is never loaded, but instead subsets of its information is queried
#   on demand utilizing the @see "PlatformDesign".
# }
# 
# \examples{\dontrun{
# }}
#
# \author{
#   Henrik Bengtsson.  
#   All credits to Benilton Carvalho and Rafael Irizarry for the implementation
#   of the underlying calibration algorithm.
# }
#
# \references{
#   Carvalho B, Bengtsson H, Speed TP and Irizarry RA. \emph{Exploration, 
#   Normalization, and Genotype Calls of High Density Oligonucleotide SNP 
#   Array Data}, 2006.\cr
#
#   Naef F and Magnasco MO. \emph{Solving the riddle of the bright 
#   mismatches: labeling and effective binding in oligonucleotide arrays}. 
#   Phys Rev E Stat Nonlin Soft Matter Phys, 2003, 68.\cr
# }
#*/###########################################################################
setConstructorS3("PcrProbeCalibration", function(..., subsetToFit=1/5, flavor=c("*", "oligo")) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  extend(ProbePreprocessing(...), "PcrProbeCalibration",
    .subsetToFit = subsetToFit,
    .flavor = flavor,
    "cached:.coefs" = NULL
  )
})

setMethodS3("getParameters", "PcrProbeCalibration", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  # Get parameters of this class
  params2 <- list(
    subsetToFit = getSubsetToFit(this)
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, private=TRUE)


setMethodS3("getPmDesignMatrix", "PcrProbeCalibration", function(this, ..., verbose=FALSE) {
  require(oligo) || throw("Package not loaded: oligo"); # gcrma_getSeq2()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting the design matrix");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the platform-design object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  pd <- PlatformDesign(cdf);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Gets the nucleotide sequences of all PM probes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting nucleotide sequences for all PMs");
  pmIdxs <- which(getFeatureInfo(pd, "feature_type") == "PM");
  nbrOfPms <- length(pmIdxs);
  # Get the number of nucleotides in each probe (typically 25)
  seqs <- getFeatureInfo(pd, "sequence", subset=pmIdxs);
  verbose && printf(verbose, "RAM: %.2fMB\n", object.size(seqs)/1024^2);
  oligoLength <- nchar(seqs[1]);
  # By, pasting 'seqs' already hear we should save at least of the RAM.
  seqs <- paste(seqs, collapse="");
  verbose && printf(verbose, "RAM: %.2fMB\n", object.size(seqs)/1024^2);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the unit names for each of the PM probe
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  unitNames <- getFeatureInfo(pd, field="feature_set_name", subset=pmIdxs);
  unitNames <- as.character(unitNames);
  rm(pmIdxs);
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Loads SNP annotations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Loading the SNP annotations");
  annot <- getAnnotations(pd);
  snpNames <- annot$SNP;
  fragLengths <- as.integer(annot$Length);
  rm(annot);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the overall design matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Building the design matrix");

  verbose && enter(verbose, "Allocating design matrix");
  size <- nbrOfPms*(1+oligoLength+3);
  verbose && printf(verbose, "Estimated RAM: %.2fMB\n", (8*size)/1024^2);
  D <- matrix(1, nrow=nbrOfPms, ncol=1+(3*oligoLength)+3);
  verbose && exit(verbose);

  verbose && enter(verbose, "Getting design-matrix elements for PCR fragment length effects");
  # Gets the PCR fragment lengths for each PM probe
  pos <- match(unitNames, snpNames);
  rm(unitNames);
  fragLengths <- fragLengths[pos];  # Expands to probe level.
  rm(pos);
  ok <- is.finite(fragLengths);
  fragLengths[!ok] <- median(fragLengths[ok]);
  rm(ok);
  verbose && print(verbose, summary(fragLengths));
  # Gets the basis matrix for natural cubic splines
  cc <- 1+3*oligoLength+1:3;
  D[,cc] <- ns(fragLengths, df=3);
  rm(fragLengths);
  verbose && exit(verbose);

  verbose && enter(verbose, "Getting design-matrix elements for sequence effects");
  # Gets the Nx(25x(4-1)) design matrix where N = #PM probes
  cc <- 1+1:(3*oligoLength);
  # Adopted from oligo::sequenceDesignMatrix().
  D[,cc] <- .Call("gcrma_getSeq2", seqs, nbrOfPms, oligoLength, PACKAGE="oligo");
  rm(seqs);
  verbose && exit(verbose);


  verbose && str(verbose, D);
  verbose && printf(verbose, "RAM: %.2fMB\n", object.size(D)/1024^2);

  verbose && exit(verbose);

  verbose && exit(verbose);

  D;
}, private=TRUE)




setMethodS3("getSubsetToFit", "PcrProbeCalibration", function(this, ...) {
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  pd <- PlatformDesign(cdf);
  pmIdxs <- which(getFeatureInfo(pd, "feature_type") == "PM");
  N <- length(pmIdxs);

  flavor <- this$.flavor;
  if (flavor == "oligo") {
    # From 'oligo': 
    set.seed(1); 
    subset <- sample(1:N, size=10000);
  } else {
    # Get the subset of PM probes to fit
    subsetToFit <- this$.subsetToFit;
    subset <- round(seq(from=1, to=N, length.out=subsetToFit*N));
  }

  subset;
}, private=TRUE)




setMethodS3("fit", "PcrProbeCalibration", function(this, D, ly, force=FALSE, ...) {
  # Adopted from oligo::CalibrateSequenceLength()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Check cache
  coefs <- this$.coefs;
  if (!force && !is.null(coefs))
    return(coefs);

  # Precalculation before model fitting to save memory
  tD <- t(D);
  D <- solve(tD %*% D);
  D <- D %*% tD;
  rm(tD);

  # Fit the model to the log2 signals
  coefs <- D %*% ly;

  colnames(coefs) <- getNames(getInputDataSet(this));
  oligoLength <- (nrow(coefs) - (1+3))/3;
  probes <- rep(c("A","C","T"), each=oligoLength);
  probes <- sprintf("%s.%02d", probes, seq(length=oligoLength));
  rownames(coefs) <- c("offset", probes, 1:3);

  # Update cache
  this$.coefs <- coefs;

  coefs;
}, private=TRUE)



setMethodS3("process", "PcrProbeCalibration", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Correcting data set for sequence effects");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already processed");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  pd <- PlatformDesign(cdf);

  # Get algorithm parameters
  params <- getParameters(this);

  # Get (and create) the output path
  outputPath <- getPath(this);

  verbose && cat(verbose, "Output path: ", outputPath);

  nbrOfArrays <- nbrOfArrays(ds);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the PMs according to oligo
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pmIdxs <- which(getFeatureInfo(pd, "feature_type") == "PM");

  # Remap oligo indices as CDF indices
  map <- getWriteMap(pd);
  pmIdxs <- map[pmIdxs];
  rm(map);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Get the design matrix for all PM probes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # This requires a lot of memory for large CDFs!
  D <- getPmDesignMatrix(this, verbose=less(verbose)); 

  # Get the subset to fit
  subset <- getSubsetToFit(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit the model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Fitting the model to a subset of the data");
  # Retrieve the the subset of data from *all* arrays to be fitted
  verbose && enter(verbose, "Loading subset of PM signals across arrays to fit the model");
  verbose && cat(verbose, "Number of arrays: ", nbrOfArrays);
  nSubset <- length(subset);
  verbose && printf(verbose, "Number of PM probes used/array: %d (%.1f%%)\n", 
                                      nSubset, 100*nSubset/length(pmIdxs));
  verbose && printf(verbose, "Required amount of RAM: > %.2fMB\n", 
                                           (8*nbrOfArrays*nSubset)/1024^2);
  pathnames <- getPathnames(ds);
  y <- readCelIntensities(pathnames, indices=pmIdxs[subset]);
  verbose && str(verbose, y);
  verbose && cat(verbose, "Transforming to log2 scale");
  # Work with the log2 signals
  y <- log2(y);
  verbose && str(verbose, y);
  verbose && exit(verbose);
  verbose && enter(verbose, "Fitting");
  coefs <- fit(this, D=D[subset,], ly=y, force=force, verbose=less(verbose));
  rm(y);
  verbose && str(verbose, coefs);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibrated each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Pre-calculations to minimize memory usage
  D <- (D %*% coefs);
  rm(coefs);

  for (kk in seq(length=nbrOfArrays)) {
    # Grab the source CEL file
    df <- getFile(ds, kk);

    verbose && enter(verbose, sprintf("Calibrating array %d of %d", kk, nbrOfArrays));
    verbose && cat(verbose, "Full name: ", getFullName(df));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Load data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Loading all PM intensities");
    y <- getData(df, indices=pmIdxs, fields="intensities")$intensities;
    verbose && str(verbose, y);
    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calibrate data on the log scale
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Calibrating on the log2 scale");
    y <- log2(y);
    # Get the average signal for each array
    # Correct for sequence effects
    y <- (y - D[,kk]) + mean(y, na.rm=TRUE);
    # Go back to the intensity scale
    y <- 2^y;
    verbose && str(verbose, y);
    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Storing the calibrated data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing data");
    filename <- sprintf("%s.cel", getFullName(df));
    pathname <- filePath(outputPath, filename);
    verbose && cat(verbose, "Pathname: ", pathname);
    # Copy CEL file and update the copy
    verbose && enter(verbose, "Copying source CEL file");
    copyCel(from=getPathname(df), to=pathname, overwrite=force);
    verbose && exit(verbose);
    # Updating CEL
    verbose && enter(verbose, "Writing calibrated PM intensities");
    updateCel(pathname, indices=pmIdxs, intensities=y);
    rm(y);
    verbose && exit(verbose);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (kk ...)  

  rm(D);

  # Update the output data set
  outputDataSet <- getOutputDataSet(this);
  setCdf(outputDataSet, cdf);
  this$outputDataSet <- outputDataSet;

  verbose && exit(verbose);
  
  outputDataSet;
})



setMethodS3("plotFit", "PcrProbeCalibration", function(this, ylim=c(-1,1), xlab="Nucleotide position", ylab="coefficient", ...) {
  # Check cache
  coefs <- this$.coefs;
  if (is.null(coefs))
    throw("Model has to be fitted first.");

  # Keep only probe parameters
  excl <- c(1,seq(from=nrow(coefs)-2,to=nrow(coefs)));
  coefs <- coefs[-excl,];
  coefs <- t(coefs);
  coefs <- unwrap(coefs); # R.utils

  cols <- c(A="red", C="blue", T="orange", G="purple");
  dimnames <- dimnames(coefs);
  oligoLength <- length(dimnames[[3]]);
  xlim <- c(1,oligoLength);
  plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab);
  legend("topright", legend=names(cols[-4]), col=cols[-4], pch=19);
  for (aa in dimnames[[1]]) {
    for (nn in dimnames[[2]]) {
      params <- coefs[aa,nn,];
      x <- seq(along=params);
      lines(x, params, col=cols[nn]);
    }
  }
  invisible(coefs);
}, private=TRUE)


############################################################################
# HISTORY:
# 2006-12-12
# o Added argument 'flavor' to constructor.
# 2006-12-08
# o Now this pre-processor output results to probeData/.
# 2006-12-07
# o Renamed to PcrProbeCalibration to indicate that it is a method 
#   correcting for PCR effects at the probe level.
# o Verified that this class gives the same corrected PM values as the
#   oligo package gives.
# o Verified that fit() returns the same as oligo would do.
# o Verified that getPmDesignMatrix() returns the same as cbind(1,SeqMat,L)
#   in oligo::preProcess().
# o Made more memory efficient by only loading subset of data to be fitted,
#   and then calibrating one array at the time.  Thus, the memory 
#   limitation is now in the fit of the model, not the calibration.  When
#   fitting the model only PM signals are used, which is roughly half the
#   number of probes on an array.  In turn, only a subset (20%) of these 
#   probes are used to fit the model, so effectively we have to have memory
#   enough to load 10% of probe signals.
# o First version running using mostly oligo calls.
# o Created from QuantileNormalizer.R.
############################################################################
