###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod isSnpChip
#
# @title "Static method to check if a chip is a mapping (SNP) chip"
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
#   Returns @TRUE if the chip type refers to a SNP array, otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isSnpChip", "AffymetrixCdfFile", function(this, ...) {
  chipType <- getChipType(this);

  if (regexpr("^Mapping(10K|50K|250K)_.*$", chipType) != -1)
    return(TRUE);

  if (regexpr("^Cent(Hind|Xba).*$", chipType) != -1)
    return(TRUE);

  FALSE;
}, private=TRUE)



###########################################################################/**
# @RdocMethod getSnpNames
#
# @title "Gets the names of the SNPs"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{pattern}{A regular expression to identify unit that are SNPs.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   Internally, @seemethod "getUnitNames".
#   is used.
#
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getSnpNames", "AffymetrixCdfFile", function(this, pattern="^SNP_A-", ...) {
  unitNames <- getUnitNames(this, ...);
  grep(pattern=pattern, unitNames, value=TRUE);
}, private=TRUE)




###########################################################################/**
# @RdocMethod nbrOfSnps
#
# @title "Gets the number of SNPs"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional arguments passed to @seemethod "getSnpNames".}
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   Internally, @seemethod "getSnpNames" is used to identify SNPs.
#
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfSnps", "AffymetrixCdfFile", function(this, ...) {
  length(getSnpNames(this, ...));
}, private=TRUE)




###########################################################################/**
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


###########################################################################/**
# @RdocMethod getAlleleProbePairs
#
# @title "Gets the indices of probepairs with the same pair of SNP nucleotides"
#
# \description{
#   @get "title".
#   Note that the order of allele A and allele B is irrelevant.
#   For instance, all probepairs with nucleotides (A,T) are calibrated
#   together with all probepairs with nucleotides (T,A) reversed.
# }
#
# @synopsis
#
# \arguments{
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @list where each element is a two-column @matrix where
#   the column names are the nucleotides for the two alleles.
# }
#
# \section{Benchmarking}{
#   On an IBM Thinkpad A31 with 1.8GHz and 1GB RAM:
#   \itemize{
#    \item{Mapping10K\_Xba142}{10208 units & 432964 cells: 11 seconds.}
#    \item{Mapping50K\_Xba240}{58960 SNPs & 589,600 (PMA,PMB) probe pairs: 11 seconds.}
#   }
# }
#
# @author
#
# \seealso{
#   @seemethod "getAllelePairUnitSets".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAlleleProbePairs", "AffymetrixCdfFile", function(this, units=NULL, ignoreOrder=TRUE, force=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Identifying the probes stratified by allele basepairs");
  on.exit(verbose && exit(verbose));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this);
  key <- list(method="getAlleleProbePairs", class=class(this)[1], chipType=chipType, units=units, ignoreOrder=ignoreOrder);
  dirs <- c("aroma.affymetrix", chipType);
  if (!force) {
    probeSets <- loadCache(key=key, dirs=dirs);
    if (!is.null(probeSets))
      return(probeSets);
  }

  cdfFile <- getPathname(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify all possible allele pairs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Loading all possible allele basepairs");
  # Use only units that are SNPs
  unitNames <- readCdfUnitNames(cdfFile);
  unitsAll <- which(regexpr("^SNP", unitNames) != -1);
  rm(unitNames);
  gc();

  # Operate only on a subset of probes?
  if (!is.null(units)) {
    unitsAll <- intersect(unitsAll, units);
  }
  units <- unitsAll;
  rm(unitsAll);
  nunits <- length(units);

  verbose && cat(verbose, "Number of units to query: ", nunits);
  if (nunits == 0)
    return(NULL);

  # Read group names for these SNPs
  verbose && enter(verbose, "Requiring group names");
  groupNames <- readCdfGroupNames(cdfFile, units=units);
  # Save memory by removing names. [55Mb -> 44Mb]
  names(groupNames) <- NULL;
  # Save memory by converting to integers. [44Mb -> 11Mb]
  levels <- as.integer(1:4);
  names(levels) <- c("A", "C", "G", "T");
  groupNames <- base::lapply(groupNames, FUN=function(s) { 
    s <- levels[s];
    names(s) <- NULL;
    s;
  });
  uGroupNames <- unique(groupNames);
  gc();
  verbose && exit(verbose);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read all of the CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Loading cell indices for probepairs for requested units");
  cdfAll <- readCdfCellIndices(cdfFile, units=units, stratifyBy="pm");
  # Save memory by removing names. [309Mb -> 298Mb]
  names(cdfAll) <- NULL;

  cells0 <- unlist(cdfAll, use.names=FALSE);
  nbrOfCells <- length(cells0);
  verbose && printf(verbose, "Identified %d (PM_A,PM_B) pairs in %d units\n", 
                                                round(nbrOfCells/2), nunits);
  cells0 <- sort(cells0);

  # Not needed anymore
  rm(units);
  gc();

  # Save memory by flattening structure. [298Mb -> 51Mb(!)]
  # TODO: Add support to do this already in affxparser?! /HB 2006-07-22
  cdfAll <- base::lapply(cdfAll, FUN=function(unit) {
    groups <- .subset2(unit, 1);
    names(groups) <- NULL;
    base::lapply(groups, FUN=.subset2, 1);
  });
  gc();
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Group all units with the same allele basepairs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Stratifying by unique allele basepairs");
  probeSets <- vector("list", length(uGroupNames));
  for (kk in 1:length(uGroupNames)) {
    name <- uGroupNames[[kk]];
    basepair <- paste(names(levels)[name[1:2]], collapse="");
    verbose && enter(verbose, "Allele basepair ", basepair);

    idx <- base::lapply(groupNames, FUN=identical, name);
    idx <- which(unlist(idx, use.names=FALSE));
    cdf <- cdfAll[idx];
    cdfAll[idx] <- NA;  # Not needed anymore
    rm(idx);

    cdf0 <- vector("list", length=length(name));
    for (gg in 1:length(name)) {
      cdf0[[gg]] <- unlist(base::lapply(cdf, FUN=.subset2, gg), use.names=FALSE);
    }
    rm(cdf);
    probeSets[[kk]] <- cdf0;
    rm(cdf0);
    names(probeSets)[kk] <- basepair;

    verbose && exit(verbose);
  }
  rm(cdfAll);
  verbose && exit(verbose);

  # Assert correctness
  verbose && enter(verbose, "Asserting correctness part I", level=-20);
  nbrOfCells2 <- length(unlist(probeSets, use.names=FALSE));
  if (nbrOfCells2 != nbrOfCells) {
    throw("Internal error1: Excepted ", nbrOfCells, " indices: ", nbrOfCells2);
  }
  if (!identical(sort(unlist(probeSets, use.names=FALSE)), cells0)) {
    throw("Internal error1: Mismatching probes.");
  }
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify equivalent groups
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Putting equivalent groups together");
  probeSets2 <- list();
  for (kk in 1:length(probeSets)) {
    bp <- names(probeSets)[kk];
    value <- probeSets[[kk]];
    for (ll in 0:(1*(length(value) > 2))) {
      value2 <- probeSets2[[bp]];
      if (is.null(value2))
        value2 <- vector("list", length=2);
      value2[[1]] <- c(value2[[1]], value[[1]]);
      value2[[2]] <- c(value2[[2]], value[[2]]);
      probeSets2[[bp]] <- value2;
      bp <- strsplit(bp, split="")[[1]];
      bp <- c(A="T", C="G", G="C", T="A")[bp];
      bp <- paste(bp, collapse="");
      value <- value[3:4];
    }
  }
  verbose && exit(verbose);

  # Assert correctness
  verbose && enter(verbose, "Asserting correctness part II", level=-20);
  nbrOfCells2 <- length(unlist(probeSets, use.names=FALSE));
  if (nbrOfCells2 != nbrOfCells) {
    throw("Internal error2: Excepted ", nbrOfCells, " indices: ", nbrOfCells2);
  }
  if (!identical(sort(unlist(probeSets, use.names=FALSE)), cells0)) {
    throw("Internal error1: Mismatching probes.");
  }
  verbose && exit(verbose);

  if (ignoreOrder) {
    verbose && enter(verbose, "Putting AB and BA groups together");
    pairs <- strsplit(names(probeSets2), split="");
    pairs <- base::lapply(pairs, FUN=function(x) paste(sort(x), collapse=""));
    pairs <- unlist(pairs);
    uPairs <- sort(unique(pairs));
    probeSets <- list();
    for (pair in uPairs) {
      idx <- which(pairs == pair);
      basepairs <- sort(names(probeSets2)[idx]);
      probeSets[[pair]] <- probeSets2[basepairs];
    }
    rm(probeSets2);
    verbose && exit(verbose);
  
    verbose && enter(verbose, "Combining AB and BA groups");
    # Join AB with BA.
    for (kk in 1:length(probeSets)) {
      values <- probeSets[[kk]];
      if (length(values) > 1) {
        values[[1]][[1]] <- c(values[[1]][[1]], values[[2]][[2]]);
        values[[1]][[2]] <- c(values[[1]][[2]], values[[2]][[1]]);
        values <- values[[1]];
        probeSets[[kk]] <- values;
      }
    }
    rm(values);
    verbose && exit(verbose);
  } else {
    probeSets <- probeSets2;
    rm(probeSets2);
  }

  # Assert correctness
  verbose && enter(verbose, "Asserting correctness part III", level=-20);
  nbrOfCells2 <- length(unlist(probeSets, use.names=FALSE));
  if (nbrOfCells2 != nbrOfCells) {
    throw("Internal error3: Excepted ", nbrOfCells, " indices: ", nbrOfCells2);
  }
  if (!identical(sort(unlist(probeSets, use.names=FALSE)), cells0)) {
    throw("Internal error1: Mismatching probes.");
  }
  verbose && exit(verbose);

  verbose && enter(verbose, "Reformatting to matrices");
  # Order indices by allele A (just for beauty)
  for (kk in 1:length(probeSets)) {
    values <- probeSets[[kk]];
    values <- matrix(c(values[[1]], values[[2]]), ncol=2);
    colnames(values) <- strsplit(names(probeSets)[kk], split="")[[1]];
    o <- order(values[,1]);
    values <- values[o,];
    probeSets[[kk]] <- values;
  }
  rm(values);
  if (isVisible(verbose, level=-20))
    verbose && str(verbose, probeSets, level=-20);
  verbose && exit(verbose);

  # Assert correctness
  verbose && enter(verbose, "Asserting correctness part IV", level=-20);
  nbrOfCells2 <- length(unlist(probeSets, use.names=FALSE));
  if (nbrOfCells2 != nbrOfCells) {
    throw("Internal error4: Excepted ", nbrOfCells, " indices: ", nbrOfCells2);
  }
  if (!identical(sort(unlist(probeSets, use.names=FALSE)), cells0)) {
    throw("Internal error: The identified set of indices for various allele probe pairs does not match the original set of cell indices.");
  }
  verbose && exit(verbose);

  # Save cache to file
  comment <- key[c("method", "class", "chipType")];
  comment <- paste(names(comment), comment, sep="=");
  comment <- paste(comment, collapse=", ");
  saveCache(probeSets, key=key, comment=comment, dirs=dirs);

  probeSets;
}, private=TRUE) # getAlleleProbePairs()




setMethodS3("getAlleleProbePairs2", "AffymetrixCdfFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "Identifying the probes stratified by allele basepairs");
  cdfFile <- getPathname(this);

  # Identify all possible allele pairs
  verbose && enter(verbose, "Loading all possible allele basepairs");
  # Use only units that are SNPs
  unitNames <- readCdfUnitNames(cdfFile);
  units <- which(regexpr("^SNP", unitNames) != -1);
  rm(unitNames);

  # Read group names for the SNPs
  groupNames <- readCdfGroupNames(cdfFile, units=units);
  uGroupNames <- unique(groupNames);
  verbose && exit(verbose);

  uGroupNames0 <- base::lapply(uGroupNames, FUN=function(x) {
    x <- matrix(x, nrow=2)
    if (ncol(x) == 2) {
      # Take the complement bases for the reverse strand
      x[,2] <- c(A="T", C="G", G="C", T="A")[x[,2]];
    }
    x;

  })

  uBasepairs0 <- base::lapply(uGroupNames0, FUN=function(x) {
    base::apply(x, MARGIN=2, FUN=sort);
  })

  uBasepairs1 <- base::lapply(uBasepairs0, FUN=function(x) {
    base::apply(x, MARGIN=2, FUN=paste, collapse="");
  })

  # Get all unique allele basepairs
  uBasepairs <- sort(unique(unlist(uBasepairs1)));
  uBasepairs <- strsplit(uBasepairs, split="");

  # Create basepairs to group names map
  map <- vector("list", length(uBasepairs));
  names(map) <- uBasepairs;
  bpIdx <- vector("list", length(uBasepairs0));
  for (bp in uBasepairs) {
    for (kk in 1:length(bpIdx)) {
      set <- uBasepairs0[[kk]];
      bpIdx[[kk]] <- kk + which(bp == set)/10;
    }
    map[[bp]] <- unlist(bpIdx);
  }

  # Read all of the CDF file
  verbose && enter(verbose, "Loading cell indices for all probepairs");
  cdfAll <- readCdfCellIndices(cdfFile, units=units, stratifyBy="pm");
  rm(units);
  verbose && exit(verbose);

  verbose && enter(verbose, "Stratifying by unique allele basepairs");
  probes <- vector("list", length(map));
  for (kk in 1:length(map)) {
    basepair <- names(map)[kk];
    verbose && enter(verbose, "Allele basepair ", basepair);

    bpIdx <- map[[kk]];
    gIdx <- as.integer(bpIdx);
    sIdx <- round(10*(bpIdx - gIdx));
    
    verbose && cat(verbose, "Located in ", length(unique(gIdx)), " group(s).");
    
    idx <- base::lapply(groupNames, FUN=identical, basepair);
    idx <- which(unlist(idx, use.names=FALSE));
    cdf <- cdfAll[idx];

    cdf0 <- vector("list", length=4);
    for (gg in 1:4) {
      cells <- applyCdfGroups(cdf, cdfGetGroups, gg);
      cells <- unlist(cells, use.names=FALSE);
      cdf0[[gg]] <- cells;
    }
    probes[[kk]] <- cdf0;
    names(probes)[kk] <- basepair;

    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  probes;
}, private=TRUE) # getAlleleProbePairs2()



############################################################################
# HISTORY:
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
