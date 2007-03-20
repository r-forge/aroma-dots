###########################################################################/**
# @set "class=ChipEffectSet"
# @RdocMethod calculateBaseline
#
# @title "Estimates the baseline signal chromosome by chromosome"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{chromosomes}{An @integer @vector specifying for which chromsosomes
#     the baseline should be estimated.  
#     If @NULL, all chromosomes are considered.}
#   \item{ploidy}{An @integer specifying the ploidy that the baseline
#     should have.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("calculateBaseline", "ChipEffectSet", function(this, chromosomes=NULL, ploidy=2, ..., force=FALSE, verbose=FALSE) {
  cdf <- getCdf(this);
  gi <- getGenomeInformation(cdf);
  allChromosomes <- getChromosomes(gi);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosomes':
  if (is.null(chromosomes)) {
    chromosomes <- allChromosomes;
  } else {
    chromosomes <- Arguments$getChromosomes(chromosomes, 
                                                range=range(allChromosomes));
  }

  # Argument 'ploidy':
  ploidy <- Arguments$getInteger(ploidy, range=c(1,8));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Estimating the baseline signals for each chromosome");

  n <- nbrOfArrays(this);
  for (chromosome in chromosomes) {
    verbose && enter(verbose, "Chromosome ", chromosome);

    verbose && enter(verbose, "Extracting the two sets of samples");
    isBaseline <- (sapply(this, getPloidy, chromosome) == ploidy);
    nB <- sum(isBaseline);
    verbose && printf(verbose, "Number of samples with ploidy %d: %d\n",
                                                              ploidy, nB);
    if (nB == 0) {
      throw("Cannot estimate baseline signals. No samples with ploidy ", ploidy, " available.");
    }

    nM <- n - nB;
    
    # Baseline samples
    csB <- extract(this, which( isBaseline));
    verbose && printf(verbose, "Baseline samples (with ploidy %d):\n", ploidy);
    verbose && print(verbose, csB);

    if (nM > 0) {
      # Mixed samples
      csM <- extract(this, which(!isBaseline));
      verbose && cat(verbose, "All other samples:");
      verbose && print(verbose, csM);
      verbose && exit(verbose);
    }
    
    verbose && enter(verbose, "Identifying units on chromosome");
    units <- getUnitsOnChromosome(gi, chromosome=chromosome);
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);
    verbose && exit(verbose);

   verbose && enter(verbose, "Identifying ucells for these units");
    cells <- getCellIndices(this, units=units);
    rm(units);
    cells <- unlist(cells, use.names=FALSE);
    cells <- sort(cells);
    verbose && cat(verbose, "Cells:");
    verbose && str(verbose, cells);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating average of baseline samples");
    csBavg <- getAverageFile(csB, indices=cells, force=force, verbose=less(verbose));
    verbose && exit(verbose);

    if (nM > 0) {
      verbose && enter(verbose, "Calculating average of all other samples");
      csMavg <- getAverageFile(csM, indices=cells, force=force, verbose=less(verbose));
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Calculating differences of the means");
    muBs <- getData(csBavg, indices=cells, fields="intensities", verbose=less(verbose))$intensities;
    rm(csBavg);
    verbose && str(verbose, muBs);
    
    verbose && cat(verbose, "Summary of log2(mu2s)");
    verbose && print(verbose, summary(log2(muBs)));
      
    if (nM > 0) {
      muMs <- getData(csMavg, indices=cells, fields="intensities", verbose=less(verbose))$intensities;
      rm(csMavg);
      verbose && str(verbose, muMs);
      
      # On the log scale, it would have been a difference
      cs <- (muBs / muMs);
      verbose && str(verbose, cs);
      verbose && cat(verbose, "Summary of log2(cs*)");
      verbose && print(verbose, summary(log2(cs)));
      
      c <- median(cs, na.rm=TRUE);
      rm(cs);
      verbose && printf(verbose, "log2(c*) = %.3f\n", log2(c));
      verbose && exit(verbose);
      
      verbose && cat(verbose, "Summary of log2(mu1s)");
      muBs2 <- muMs * c;
      rm(muMs);
      verbose && print(verbose, summary(log2(muBs2)));

      wB <- nB/n;
      wM <- 1-wB;

      verbose && cat(verbose, "Summary of log2(ds)");
      ds <- wB*muBs + wM*muBs2;
      rm(muBs, muBs2);
    } else {
      ds <- muBs;
    }
    verbose && print(verbose, summary(log2(ds)));

    verbose && enter(verbose, "Updating baseline signals");
    # TO DO
    rm(ds);
    verbose && exit(verbose);

    rm(cells);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
    
    verbose && exit(verbose);
  }


  verbose && exit(verbose);

  res;
}) # calculateBaseline()


############################################################################
# HISTORY:
# 2007-03-16
# o Created.  See ploidy4.ps paper.
############################################################################
