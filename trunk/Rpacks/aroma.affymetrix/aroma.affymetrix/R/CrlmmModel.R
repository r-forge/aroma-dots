###########################################################################/**
# @RdocClass CrlmmModel
#
# @title "The CrlmmModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the CRLMM model of Carvalho et al. (2006).
# }
# 
# @synopsis
#
# \arguments{
#   \item{ces}{A @see "SnpChipEffectSet" to which this model should 
#     be fitted.}
#   \item{minLLRforCalls}{A @numeric @vector of length three.}
#   \item{transform}{The transform used for the chip effects.}
#   \item{tags}{A @character @vector of tags to be appended to the tags of
#      the input data set.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Requirments}{
#   
#   This implementation is partly adopted from and partly makes use of
#   the implementation available in the \pkg{oligo} package.
#   It was implemented so it can give identical results to the \pkg{oligo}
#   implementation.
# }
#
# @author
#
# \references{
#  Carvalho B, Bengtsson H, Speed TP, and Irizarry RA. \emph{Exploration, 
#  Normalization, and Genotype Calls of High Density Oligonucleotide SNP 
#  Array Data}, Biostatistics, 2006.\cr
# }
#*/###########################################################################
setConstructorS3("CrlmmModel", function(ces=NULL, minLLRforCalls=c(AA=50, AB=40, BB=50), transform=c("log", "asinh"), tags="CRLMM", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ces':
  if (!is.null(ces)) {
    if (!inherits(ces, "SnpChipEffectSet"))
      throw("Argument 'ces' is not an SnpChipEffectSet object: ", class(ces));
  }

  # Arguments 'minLLRforCalls':
  minLLRforCalls <- Arguments$getDoubles(minLLRforCalls, length=3);

  # Argument 'transform':
  transform <- match.arg(transform);

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
  }


  extend(Model(dataSet=ces, tags=tags), "CrlmmModel",
    minLLRforCalls = minLLRforCalls,
    transform = transform
  )
})



setMethodS3("as.character", "CrlmmModel", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("SNP chip-effect set: %s", getName(this)));
  ces <- getChipEffects(this);
  tags <- paste(getTags(ces), collapse=",");
  s <- c(s, sprintf("Input tags: %s", tags));
  s <- c(s, sprintf("Transform: %s", this$transform));
  s <- c(s, sprintf("Output tags: %s", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(this))));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


###########################################################################/**
# @RdocMethod getRootPath
#
# @title "Gets the root path of this model"
#
# \description{
#  @get "title", which is the class name of the model preceeded by the
#  string \code{"model"}.
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
setMethodS3("getRootPath", "CrlmmModel", function(this, ...) {
  "genotypeData";
}, private=TRUE)




###########################################################################/**
# @RdocMethod getChipEffects
#
# @title "Gets the chip effects for this model"
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
#   Returns an @see "SnpChipEffectSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChipEffects", "CrlmmModel", function(this, ...) {
  getDataSet(this);
})



###########################################################################/**
# @RdocMethod classifyAsXorXX
#
# @title "Classifies the samples to have one or two copies of chromosome X"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @vector of length N, where N is the number of samples.
# }
#
# \details{
#  This method is adopted from the internal \code{snpGenderCall()} function
#  in the \pkg{oligo} package.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("classifyAsXorXX", "CrlmmModel", function(this, ...) {
  # Setup
  ds <- getInputDataSet(this);
  cdf <- getCdf(ds);
  gi <- getGenomeInformation(cdf);

  # Get the median SNP signal across all chromosomes and arrays
  A <- getA(object);  # (SNP, sample, strand)
  medA <- median(A, na.rm=TRUE);

  # Get the median chromosome X signal for each array
  xUnits <- getUnitsOnChromosome(gi, "X");
  xUnits <- getChrXIndex(object);
  A <- A[xUnits,,,drop=FALSE];
  a <- apply(A, MARGIN=2, FUN=median, na.rm=TRUE);

  # Cluster by K-mean algorithm
  minA <- min(a, na.rm=TRUE);
  fit <- kmeans(a, centers=c(minA, medA));

  # Identify class
  class <- as.integer(kfit$cluster==1) + 1;
  class <- factor(c("XX","X")[class]);

  # Return
  class;
}, private=TRUE)


###########################################################################/**
# @RdocMethod fit
#
# @title "Fits the model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "oligo::SnpCallSetPlus-class" object.
# }
#
# @author
#
# \seealso{
#   @see "oligo::crlmm".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("fit", "CrlmmModel", function(this, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Requirements
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(oligo) || throw("Package 'oligo' not loaded.");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Default settings
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  recalibrate <- TRUE;
  returnParams <- TRUE;

#  returnParams <- FALSE;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  verbose2 <- as.logical(verbose);


  verbose && enter(verbose, "Fitting CRLMM");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get output path
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- getPath(this);
  mkdirs(path);

  paramTags <- NULL;
  if (recalibrate) {
    paramTags <- c(paramTags, "recalib");
  }
  paramTags <- paste(paramTags, collapse=",");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check if already fitted
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- sprintf("%s,crlmm%s.xdr", getFullName(this), paramTags);
  pathname <- filePath(path, filename);
  if (isFile(pathname)) {
    verbose && cat(verbose, "Already fitted");
    verbose && enter(verbose, "Loading results from file");
    verbose && cat(verbose, "Pathname: ", pathname);
    fit <- loadObject(pathname);
    verbose && exit(verbose);
    verbose && exit(verbose);
    return(fit);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for platform-design package
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert that the PD enviroment package for the CDF is installed
  cdf <- getCdf(this);
  chipType <- getChipType(cdf);
  chipType <- gsub("[,-]monocell", "", chipType);
  cdfMain <- AffymetrixCdfFile$fromChipType(chipType);
  pd <- PlatformDesign(cdfMain);
  if (!isInstalled(pd)) {
    throw("Platform-design package not installed: ", getPackageName(pd));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract chip effects as a SnpQSet object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting chip effects as a SnpQSet object");
  ces <- getChipEffects(this);
  qs <- extractSnpQSet(ces, transform=this$transform, verbose=less(verbose));
  snpNames <- featureNames(qs);
  units <- indexOf(cdf, names=snpNames);  # The units according to the CDF
  verbose && print(verbose, qs);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Dimension of data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSnps <- length(snpNames);
  nbrOfSamples <- nbrOfArrays(ces);
  verbose && printf(verbose, "Number of SNPs: %d\n", nbrOfSnps);
  verbose && printf(verbose, "Number of samples: %d\n", nbrOfSamples);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify units on chromosome X
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying units on chromosome X");

  # getChrXIndex() is robust against SNP reordering /HB 2006-10-03a
  annot <- getAnnotations(pd);
  annot <- annot[match(snpNames, annot$SNP),];
  xIndex <- which(annot$Chromosome == "chrX");

  # Alternative not requiring PD environment annotations.
  gi <- getGenomeInformation(cdf);
  xUnits <- getUnitsOnChromosome(gi, "X");         # As index by the CDF.
  xUnitNames <- getUnitNames(cdf, units=xUnits);
  xIndex2 <- which(snpNames %in% xUnitNames);
  stopifnot(identical(xIndex2, xIndex));
  verbose && str(verbose, xIndex);
  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get gender
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving genders");
  # Got gender information?
  if(is.null(qs$gender)){
    verbose && enter(verbose, "Inferring gender from chromosome X signals");
    # snpGenderCall() is robust against SNP reordering /HB 2006-10-03
    qs$gender <- oligo::snpGenderCall(qs);  # Uses 'oligo'
    verbose && exit(verbose);
  }
  maleIndex <- (qs$gender == "male");
  malesStr <- paste(which(maleIndex), " (", getNames(ces)[maleIndex], ")", sep="");
  malesStr <- paste(malesStr, collapse=", ");
  verbose && cat(verbose, "Samples classified as males (=one copy of X): ", malesStr);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get M corrections
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving M corrections");
  # It *looks like* this is robust against SNP reordering /HB 2006-10-03
  filename <- sprintf("logRatioCorrections.xdr");
  pathname <- filePath(path, filename);
  if (isFile(pathname)) {
    verbose && enter(verbose, "Correction already available on file");
    verbose && cat(verbose, "Pathname: ", pathname);
    correction <- loadObject(pathname);
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "Calculating log-ratio corrections using fitAffySnpMixture()");
    correction <- oligo::fitAffySnpMixture(qs, verbose=verbose2);
    verbose && exit(verbose);
    verbose && enter(verbose, "Saving corrections to file");
    verbose && cat(verbose, "Pathname: ", pathname);
    saveObject(correction, file=pathname);
    verbose && exit(verbose);
  }
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Loading CRLMM parameters");
  verbose && cat(verbose, "Package: ", getPackageName(pd));
  env <- getCrlmmInfo(pd);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the indices of SNPs that are not called by the HapMap project 
  # or that were poorly called.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying poor (or missing) HapMap calls according to the platform-design package");
  poorCalls <- (!env$hapmapCallIndex | env$badCallIndex);
  poorCalls <- poorCalls[match(snpNames, names(env$badRegions))];
  poorCalls <- which(poorCalls);
  verbose && cat(verbose, "Among the ", nbrOfSnps, " SNPs, there are ", length(poorCalls), " SNPs that were not called by the HapMap project or were poorly called by it.");
  verbose && cat(verbose, "Feature names:");
  verbose && str(verbose, snpNames[poorCalls]);
  verbose && exit(verbose);

  # Get rid of large enviroment again
  priors <- env$priors;
  params <- env$params;
  rm(env); # Not needed anymore

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get initial calls for poor HapMap calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting initial calls for the poor HapMap calls from the corrected M values");
  calls <- matrix(NA, nrow=nbrOfSnps, ncol=nbrOfSamples);
  calls[poorCalls,] <- oligo::getInitialAffySnpCalls(correction, 
                                       subset=poorCalls, verbose=verbose2);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get SNP genotype region parameters for poor calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving initial genotype regions for poor calls");
  # Check for already recalibrated genotype regions
  filename <- sprintf("genotypeRegions,poorCalls.xdr");
  pathname <- filePath(path, filename);
  if (isFile(pathname)) {
    rm(calls); # Not needed anymore
    verbose && enter(verbose, "Regions available on file");
    verbose && cat(verbose, "Pathname: ", pathname);
    regions <- loadObject(pathname);
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "Calculating region parameters");

    verbose && enter(verbose, "Estimate regions from corrected ratios and initial calls");
    regions <- oligo::getAffySnpGenotypeRegionParams(qs, calls,
                        correction$fs, subset=poorCalls, verbose=verbose2);
    rm(calls); # Not needed anymore
    verbose && exit(verbose);

    # Shrink genotype regions toward the priors
    verbose && enter(verbose, "Shrinking regions toward the priors");
    regions <- oligo::updateAffySnpParams(regions, priors, verbose=verbose2);
    verbose && exit(verbose);

    verbose && enter(verbose, "Saving regions to file");
    verbose && cat(verbose, "Pathname: ", pathname);
    saveObject(regions, file=pathname);
    verbose && exit(verbose);

    verbose && exit(verbose);
  }
  verbose && exit(verbose);


  verbose && enter(verbose, "Update SNP parameters for poor calls using initial genotype regions");
  params <- oligo::replaceAffySnpParams(params, regions, subset=poorCalls);
  rm(poorCalls);
  verbose && exit(verbose);

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate call distances
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating call distances");
  # SNPs in qs, params, and correction$fs must be ordered the same.
  # myDist is a JxIx3x2 matrix where J=#SNPs, I=#arrays, with
  # 3=#genotypes and # 2=#strands.
  myDist <- oligo::getAffySnpDistance(qs, params, correction$fs);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the calls from the new distances
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Updating the calls from the new distances");
  # This takes time: This is worth optimizing.
  calls <- oligo::getAffySnpCalls(myDist, xIndex, maleIndex, verbose=verbose2);
  # A JxI matrix; one call per sample and SNP.
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate log-likehood ratios
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating log-likehood ratios");
  llr <- oligo::getAffySnpConfidence(myDist, calls, xIndex, maleIndex, 
                                                         verbose=verbose2);
  # A JxI matrix; one LLR per sample and SNP.
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Recalibrate?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (recalibrate) {
    rm(myDist, regions); # Not needed anymore

    verbose && enter(verbose, "Recalibrating");

    # Check for already recalibrated genotype regions
    filename <- sprintf("genotypeRegions,recalib.xdr");
    pathname <- filePath(path, filename);
    if (isFile(pathname)) {
      rm(calls, llr); # Not needed anymore
      verbose && enter(verbose, "Recalibrated genotype regions available on file");
      verbose && cat(verbose, "Pathname: ", pathname);
      regions <- loadObject(pathname);
      verbose && exit(verbose);
    } else {
      verbose && enter(verbose, "Recalibrating genotype regions");
      verbose && enter(verbose, "Resetting poor genotype calls");
      # For each genotype class
      for (kk in 1:3) {
        bad <- which(calls == kk & llr < this$minLLRforCalls[kk]);
        calls[bad] <- NA;
        rm(bad);
      }
      verbose && cat(verbose, "Number of poor calls: ", sum(is.na(calls)),
                                               " out of ", length(calls));
      rm(llr); # Not needed anymore
      verbose && exit(verbose);

      verbose && enter(verbose, "Getting SNP genotype region parameters");
      regions <- oligo::getAffySnpGenotypeRegionParams(qs, calls, 
                                          correction$fs, verbose=verbose2);
      rm(calls); # Not needed anymore
      verbose && exit(verbose);

      verbose && enter(verbose, "Updating the calls from the new distances");
      regions <- oligo::updateAffySnpParams(regions, priors);
      verbose && exit(verbose);

      verbose && enter(verbose, "Saving recalibrated genotype regions to file");
      verbose && cat(verbose, "Pathname: ", pathname);
      saveObject(regions, file=pathname);
      verbose && exit(verbose);

      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Calculating call distances");
    myDist <- oligo::getAffySnpDistance(qs, regions, correction$fs,
                                                         verbose=verbose2);
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Updating the calls from the new distances");
    calls <- oligo::getAffySnpCalls(myDist, xIndex, maleIndex, 
                                                         verbose=verbose2);
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Calculating log-likehood ratios");
    llr <- oligo::getAffySnpConfidence(myDist, calls, xIndex, maleIndex,
                                                         verbose=verbose2);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # if (recalibrate)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Saving calls, one file per array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Saving calls, LLR confidence scores, and all distances to file");
  sampleNames <- getNames(ces);

  # Turn the 'call' column into a genotype factor
  as.GenotypeFactor <- function(x, ...) {
    # Turn into a factor
    levels <- c("-"=0, AA=1, AB=2, BB=3, NC=4);
    factor(x, levels=levels, labels=names(levels), ordered=FALSE);
  }

  # Allocate an Ix3x2 array to hold all distances on one array
  dim <- c(nbrOfUnits(cdf),3,2);
  dimnames <- dimnames(myDist[,1,,]);
  dimnames <- c(list(NULL), dimnames[-1]);
  dists <- array(NA, dim=dim, dimnames=dimnames);

  # For each array...
  for (kk in seq(length=nbrOfSamples)) {
    sampleName <- sampleNames[kk];
    verbose && enter(verbose, sprintf("Array %d ('%s')", kk, sampleName));

    verbose && enter(verbose, "Calls");
    filename <- sprintf("%s,calls.xdr", sampleName);
    pathname <- filePath(path, filename);
    verbose && cat(verbose, "Pathname: ", pathname);
    data <- as.GenotypeFactor(calls[,kk]);
    attr(data, "sampleName") <- sampleName;  # Update the sample name
    saveObject(data, file=pathname);
    verbose && exit(verbose);

    verbose && enter(verbose, "Confidence scores");
    filename <- sprintf("%s,score.xdr", sampleName);
    pathname <- filePath(path, filename);
    verbose && cat(verbose, "Pathname: ", pathname);
    data <- llr[,kk];
    attr(data, "sampleName") <- sampleName;  # Update the sample name
    saveObject(data, file=pathname);
    verbose && exit(verbose);

    verbose && enter(verbose, "Distances to all genotype groups");
    filename <- sprintf("%s,distances.xdr", sampleName);
    pathname <- filePath(path, filename);
    verbose && cat(verbose, "Pathname: ", pathname);
    dists[units,,] <- myDist[,kk,,];
    attr(dists, "sampleName") <- sampleName;  # Update the sample name
    saveObject(dists, file=pathname);
    verbose && exit(verbose);

    verbose && exit(verbose);
  }
  rm(data, dists, filename, pathname, as.GenotypeFactor); # Not needed anymore
  verbose && exit(verbose);
  rm(myDist); # Not needed anymore


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculating corrected log-ratios
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Correcting log-ratios for SNP called to be homozygote");
  # Allocate for all units in the CDF
  dimnames <- list(NULL, c("antisense", "sense"));
  M <- array(NA, dim=c(nbrOfUnits(cdf),2), dimnames=dimnames);
  SIGN <- c(+1,0,-1);
  # For each sample j=1,...,J
  for (jj in seq(length=nbrOfSamples)) {
    sampleName <- sampleNames[jj];
    verbose && enter(verbose, sprintf("Array %d ('%s')", jj, sampleName));

    # The the antisense and sense log-ratios
    # Should ideally be retrieved from the 'ces' object! /HB 2006-12-13
    Ma <- assayData(qs)$antisenseThetaA[,jj]-assayData(qs)$antisenseThetaB[,jj];
    Ms <- assayData(qs)$senseThetaA[,jj]-assayData(qs)$senseThetaB[,jj];

    # For the two homozygote groups k=(1,3)=(AA,BB)
    for(kk in c(1,3)) {
      sgn <- SIGN[kk];
      # Identify SNPs with genotype kk
      iis <- which(calls[,jj] == kk);
      delta <- sgn*(regions$f0 - correction$fs[iis,jj,1]);
      Ma[iis] <- Ma[iis] + delta;
      delta <- sgn*(regions$f0 - correction$fs[iis,jj,2]);
      Ms[iis] <- Ms[iis] + delta;
    }

    verbose && enter(verbose, "Saving corrected log ratios to file");
    filename <- sprintf("%s,cM.xdr", sampleName);
    pathname <- filePath(path, filename);
    M[units,1] <- Ma;
    M[units,2] <- Ms;
    saveObject(M, file=pathname);
    verbose && exit(verbose);

    verbose && exit(verbose);
  }
  rm(SIGN, sgn, iis, delta, M); # Not needed anymore
  verbose && exit(verbose);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 
  # Anything below is actually redundant and can be retrieved/calculated
  # from the data files stored in the output directory! /HB 2006-12-13
  # 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Feature data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating feature-data object");
  if (returnParams) {
    regions <- as.data.frame(regions);
    featureData <- new("AnnotatedDataFrame",
      data=regions,
      varMetadata=data.frame(labelDescription=colnames(regions),
                             row.names=colnames(regions))
    );
  } else {
    featureData <- new("AnnotatedDataFrame");
  }
  rm(regions);  # Not needed anymore
  verbose && print(verbose, featureData);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Phenotype data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating phenotype result");
  gender <- rep("female", nbrOfSamples);
  gender[maleIndex] <- "male";
  rm(maleIndex);  # Not needed anymore
  if (is.null(qs$gender)) {
    phenoData <- new("AnnotatedDataFrame",
                        data=cbind(pData(qs),
                          data.frame(crlmmSNR=as.numeric(correction$snr),
                                     gender=gender,
                                     row.names=getNames(ces))),
                        varMetadata= rbind(varMetadata(qs),
                          data.frame(labelDescription=c("crlmmSNR", "gender"),
                                     row.names=c("crlmmSNR", "gender"))))
  } else {
    phenoData <- new("AnnotatedDataFrame",
                        data=cbind(pData(qs),
                          data.frame(crlmmSNR=as.numeric(correction$snr),
                                     row.names=getNames(ces))),
                        varMetadata= rbind(varMetadata(qs),
                          data.frame(labelDescription=c("crlmmSNR"),
                                     row.names=c("crlmmSNR"))))
  }
  rm(correction);  # Not needed anymore
  verbose && print(verbose, phenoData);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return genotype estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating SnpCallSetPlus object to return");
  fit <- new("SnpCallSetPlus",
           featureData=featureData,
           phenoData=phenoData,
           experimentData=experimentData(qs),
           annotation=annotation(qs),
           calls=calls,
           callsConfidence=llr,
           logRatioAntisense=NULL,
           logRatioSense=NULL,
         );
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Saving estimates to file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  filename <- sprintf("%s,crlmm%s.xdr", getFullName(this), paramTags);
  pathname <- filePath(path, filename);
  verbose && enter(verbose, "Saving to file");
  verbose && cat(verbose, "Pathname: ", pathname);
  saveObject(fit, file=pathname);
  verbose && exit(verbose);

  verbose && exit(verbose);

  fit;
})

############################################################################
# HISTORY:
# 2007-01-06
# o CrlmmModel now inherits from Model.
# o Renamed to CrlmmModel (from Crlmm).
# 2006-12-20
# o Added default tag "CRLMM".
# o Changed the root path from crlmmData/ to genotypeData/.
# o The filename for LLR confidence scores was changed to *,scores.xdr.
# 2006-12-16
# o Now the calls and the LLRs are stored in seperate files.
# 2006-12-13
# o TO DO: Create GenotypeCallFile and GenotypeCallSet classes, with
#   subclasses CrlmmGenotypeCallFile etc. Add extractSnpCallSetPlus()
#   to GenotypeCallSet.
# o Now the following are stored to binary (xdr) files:
#   One file per sample (expanded to one row per CDF unit):
#    a) calls, confidence scores are stored.
#    b) distances to all genotype regions for each SNP.
#    c) corrected log2(thetaA/thetaB) ratios (redundant).
#   One file per data set:
#    a) Correction parameters for the log2(A/B) ratios.
#    b) Parameters for all genotype regions for all SNPs (oligo indices)
#    c) Estimated SNP regions for SNPs with poor HapMap calls.
# o Speed up of the code calculating correction log ratios; now the code
#   is vectorizing the SNPs (and not the samples), or in other words, 
#   now the code is looping over the samples, which is much faster than
#   looping over the SNPs (which are typically of order of magnitude more).
# o It finally works, that is, fit() can now return a valid SnpCallSetPlus.
# 2006-12-12
# o Trying to get fit2() to work. Still some problem creating a valid
#   SNPCallSetPlus object at the end.  Use fit() for now.
# 2006-11-29
# o Now identification of chromosome X units can also be done using 
#   aroma.affymetrix GenotypeInformation (e.g. from dChip).  Currently both
#   the oligo and the aroma.affymetrix methods are used to validate each
#   other.
# 2006-11-20
# o Converted fitCrlmm() for SnpChipEffectSet into a model class.
# 2006-10-03
# o Created based on oligo::crlmm() v0.99.20.
############################################################################
