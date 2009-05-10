setMethodS3("extractDataForSegmentation", "RawGenomicSignals", function(this, order=TRUE, dropZeroWeights=TRUE, dropWeightsIfAllEqual=TRUE, defaultChromosome=0L, defaultSampleName="Unnamed sample", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extracting data used by segmentaion algorithms");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  sampleName <- this$fullname;
  if (is.null(sampleName)) {
    sampleName <- defaultSampleName;
  }

  chromosome <- as.integer(this$chromosome);
  if (is.na(chromosome)) {
    chromosome <- defaultChromosome;
  }
  nbrOfLoci <- nbrOfLoci(this);
  verbose && cat(verbose, "Sample name: ", sampleName);
  verbose && cat(verbose, "Chromosome: ", chromosome);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  pos <- getPositions(this);
  verbose && cat(verbose, "Positions:");
  verbose && str(verbose, pos);
  verbose && summary(verbose, pos);
  # Sanity check
  stopifnot(length(pos) == nbrOfLoci);

  signals <- getSignals(this);
  verbose && cat(verbose, "Signals:");
  verbose && str(verbose, signals);
  verbose && summary(verbose, signals);
  # Sanity check
  stopifnot(length(signals) == nbrOfLoci);

  weights <- getWeights(this);
  hasWeights <- (length(weights) > 0);
  if (hasWeights) {
    verbose && cat(verbose, "Weights:");
    verbose && str(verbose, weights);
    verbose && summary(verbose, weights);
    # Sanity check
    stopifnot(length(weights) == nbrOfLoci);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Order along genome?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (order) {
    verbose && enter(verbose, "Order data along genome");
    # Order signals by their genomic location
    o <- order(pos);
    pos <- pos[o];
    signals <- signals[o];  
    if (hasWeights) {
      weights <- weights[o];
    }
    verbose && str(verbose, pos);
    verbose && str(verbose, signals);
    if (hasWeights) {
      verbose && str(verbose, weights);
    }
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (hasWeights && dropZeroWeights) {
    # Dropping loci with non-positive weights
    verbose && enter(verbose, "Dropping loci with non-positive weights");
    keep <- whichVector(weights > 0);
    pos <- pos[keep];
    signals <- signals[keep];
    weights <- weights[keep];
    nbrOfLoci <- length(pos);
    rm(keep);
    verbose && cat(verbose, "Number of loci dropped: ", 
                                             nbrOfLoci(this)-nbrOfLoci);

    # Sanity check
    if (nbrOfLoci == 0) {
      throw("No loci with non-positive weights remains.");
    }
    verbose && exit(verbose);
  }


  if (hasWeights && dropWeightsIfAllEqual) {
    # Are all weights equal?
    verbose && enter(verbose, "Checking if all (remaining) weights are identical");
    t <- weights - weights[1];
    if (all(isZero(t))) {
      verbose && cat(verbose, "Dropping weights, because all weights are equal: ", weights[1]);
      hasWeights <- FALSE;
    }
    rm(t);
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Setting up data frame");
  data <- data.frame(chromosome=chromosome, position=pos, signal=signals);
  if (hasWeights) {
    data$weight <- weights;
  }
  attr(data, "sampleName") <- sampleName;
  verbose && str(verbose, data);
  verbose && exit(verbose);
  
  verbose && exit(verbose);

  data;
}, protected=TRUE)

############################################################################
# HISTORY:
# 2009-05-10
# o This method supports all the segmentByNnn() methods.
# o Created.
############################################################################
