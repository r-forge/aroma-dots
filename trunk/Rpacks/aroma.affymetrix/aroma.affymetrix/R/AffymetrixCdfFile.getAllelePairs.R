###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod getAllelePairs
#
# @title "Gets the allele pairs for each unit"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{units}{Subset of units to be queried.  If @NULL, all units are used.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @vector of @factors.
# }
#
# @author
#
# \seealso{
#   @seemethod "getAllelePairUnitSets".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAllelePairs", "AffymetrixCdfFile", function(this, units=NULL, ...) {
  # Generate all possible allele pairs (on the forward strand).
  # Note that the 100K and the 500K arrays only have six of these,
  # but to be sure the code will work in the future too, we define all.
  bases <- c("A", "C", "G", "T");
  pairs <- outer(bases, bases, FUN=paste, sep="");
  pairs <- matrix(pairs, nrow=4, ncol=4, byrow=TRUE);
  diag(pairs) <- NA;
  pairs <- as.vector(pairs);
  pairs <- na.omit(pairs);

  # Get the group names for all the SNPs, which are the same
  # as the allele target bases
  pathname <- getPathname(this);
  bases <- readCdfGroupNames(pathname, units=units, ...);
  
  # Build a table of allele pairs on the forward strand.  The alleles
  # on the reverse strand are complementary, if that strand exists.
  forward <- base::lapply(bases, FUN=function(n) paste(n[1:2], collapse=""));
  unitNames <- names(forward);
  forward <- unlist(forward, use.names=FALSE);
  names(forward) <- unitNames;

  factor(forward, levels=pairs);
}, private=TRUE) # getAllelePairs()



###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod getAllelePairUnitSets
#
# @title "Gets the indices of units for all possible allele pairs"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getAllelePairs".}
# }
#
# \value{
#   Returns a named @list.
# }
#
# @author
#
# \seealso{
#   @seemethod "getAllelePairs".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAllelePairUnitSets", "AffymetrixCdfFile", function(this, ...) {
  # Get the CDF file for this SNP data file
  allelePairs <- getAllelePairs(this, ...);

  levels <- levels(allelePairs);
  res <- vector("list", length=length(levels));
  names(res) <- levels;

  for (level in levels) {
    units <- which(allelePairs == level);
    if (length(units) > 0) {
      res[[level]] <- units;
    }
  } # for (pair ...)

  res;
}, private=TRUE)




############################################################################
# HISTORY:
# 2008-05-10
# o ROBUSTNESS: Added backward compatibility for cases when the cached 
#   results has sets$nonSNPs as a list.
# 2008-03-26
# o CLEAN UP: getAlleleProbePairs() of AffymetrixCdfFile would print *all*
#   identified non-SNP cells in the verbose output, instead of using str().
# 2008-02-27
# o Now getAlleleProbePairs() also returns element 'nonSNPs' (unless NULL).
# 2008-02-21
# o Now getAlleleProbePairs() only consider SNPs with 2 or 4 groups, because
#   at least one custom SNP chip we've seen a few SNPs with also 6 groups
#   (which turned out to all have the same direction).
# o GENERALIZED: Now getSnpNames(), getCnNames(), getAlleleProbePairs(),
#   getAlleleProbePairs2(), and isSnpChip() all infer unit type (SNP or CN)
#   from the CDF unit type and no longer from the unit names.
# 2007-09-14
# o Added getCnNames().
# o Updated isSnpChip() to recognize 5.0 and 6.0 chips.
# o Update regular expression for getSnpNames().
# 2007-08-16
# o Now getAlleleProbePairs() of AffymetrixCdfFile processes the CDF in
#   chunks in order to save memory.  Before the GenomeWideSNP_6 CDF would
#   consume 1.5-2.0GB RAM, but now it is using less than 500MB.
# 2007-06-11
# o BUG FIX: getAlleleProbePairs2() used non-existing object 'name' instead
#   of 'basepair'.  getAlleleProbePairs2() is currently not used anyway.
# 2006-09-15
# o Adopted to the new aroma.affymetrix structure.
# 2006-07-21
# o Added getAllelePairProbes().
# o Added getAllelePairProbesets().
# 2006-06-04
# o Added getAllelePairs().
# 2006-05-31
# o Added getSnpNames() and nbrOfSnps().
# 2006-05-30
# o Added static fromFile() which tries to call ditto of all subclasses.
# o Added static isSnpChip().
# 2006-03-30
# o Updated according to affxparser.
# 2006-03-27
# o Added detailed Rdoc comments to getRelativeAlleleSignals().
# 2006-03-24
# o Added references to DM articles and Affymetrix manuals.
# o Further speed up by improve rearrangement of CDF structure. Now a Hind
#   chip takes about 11-13 minutes instead.  11 minutes compared with 
#   35 hours is 190 times faster.
# o After several speed improvements (also in affxparser), estimation of DM
#   rank scores now takes about 15-18 minutes for the 100K Hind chip.
#   The first draft took 30-35 hours(!) and yesterday 60-80 minutes.  Note,
#   the first draft was not "stupid" code; there is always room for 
#   improvement.
# o Defined a local colSums() in getDmRankScores() specialized for matrices.
#   The overhead of the default colSums() is about 50%.
# 2006-03-23
# o Moved all SNP related methods into the new class AffymetrixSnpCelFile.
# o Added getRelativeAlleleSignals().  Note, it was designed to be used 
#   with the 10K SNP chips.  These are designed so that there are equal
#   number of forward and reverse quartets with matching offsets in both
#   strands.  This is not the case for the 100K chips and above.
# o Created.
############################################################################
