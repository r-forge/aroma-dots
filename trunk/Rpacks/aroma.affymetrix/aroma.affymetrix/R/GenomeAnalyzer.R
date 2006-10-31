###########################################################################/**
# @RdocClass GenomeAnalyzer
#
# @title "The GenomeAnalyzer class"
#
# \description{
#  @classhierarchy
#
#  This abstract class provides basic methods for any type of genome analysis.
#  Specific analysis, such as copy-number analysis or exon analysis, is
#  done by subclasses of this class.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "SnpProbeAffinityFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
#*/###########################################################################
setConstructorS3("GenomeAnalyzer", function(dataSet=NULL, ...) {
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixCelSet")) {
      throw("Argument 'dataSet' is a non-AffymetrixCelSet object: ", 
                                                              class(dataSet));
    }
  }

  extend(Object(...), "GenomeAnalyzer",
    .dataSet = dataSet
  )
})

setMethodS3("as.character", "GenomeAnalyzer", function(this, ...) {
})


setMethodS3("getDataSet", "GenomeAnalyzer", function(this, ...) {
  this$.dataSet;
})

setMethodS3("nbrOfArrays", "GenomeAnalyzer", function(this, ...) {
  ds <- getDataSet(this);
  nbrOfArrays(ds);
})


setMethodS3("getChipType", "GenomeAnalyzer", function(this, ...) {
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  getChipType(cdf);
})


###########################################################################/**
# @RdocMethod getChromosomes
#
# @title "Gets the chromosomes available"
#
# \description{
#  @get "title".# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChromosomes", "GenomeAnalyzer", function(static, ...) {
  c(1:22, "X");
}, static=TRUE)



setMethodS3("getGenomeInformation", "GenomeAnalyzer", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving genome informations");
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  gi <- getGenomeInformation(cdf);
  verbose && exit(verbose);

  gi;
})


###########################################################################/**
# @RdocMethod getModel
#
# @title "Gets the model used by this analyzer"
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
#  Returns a model.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getModel", "GenomeAnalyzer", abstract=TRUE);



###########################################################################/**
# @RdocMethod getUnitsTodo
#
# @title "Gets information about units per chromosome"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{chromosomes}{The chromosomes of interest.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a named @list.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getUnitsTodo", "GenomeAnalyzer", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  unitsTodo <- this$.unitsTodo;
  if (force || is.null(unitsTodo)) {
    model <- getModel(this);
    unitsTodo <- list();
    for (chr in getChromosomes(this)) {
      verbose && enter(verbose, "Chromosome ", chr);
  
      units <- getUnitIndices(gi, chromosome=chr);
      verbose && printf(verbose, "Found in total %d units on chromosome", length(units));
      verbose && str(verbose, units);
  
      utodo <- findUnitsTodo(model, units=units, verbose=less(verbose));
      verbose && printf(verbose, "Out of which %d units are not fitted", length(utodo));
  
      unitsTodo[[chr]] <- list(
        units = units,
        todo = utodo
      );
      verbose && exit(verbose);
    } # for (chr in ...)

    this$.unitsTodo <- unitsTodo;
  }

  unitsTodo;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getProgress
#
# @title "Gets the fraction of units process per chromosome"
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
#  Returns a named @list.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getProgress", "GenomeAnalyzer", function(this, ...) {
  unitsTodo <- getUnitsTodo(this, ...);
  todo <- unlist(lapply(unitsTodo, FUN=function(chr) {
    length(chr$todo);
  }))
  counts <- unlist(lapply(unitsTodo, FUN=function(chr) {
    length(chr$units);
  }))
  df <- data.frame(counts=counts, todo=todo, progress=(counts-todo)/counts);
  total <- colSums(df);
  total["progress"] <- 1-total["todo"]/total["counts"]
  df <- rbind(df, total);
  rownames(df)[nrow(df)] <- "total";
  df;
})



###########################################################################/**
# @RdocMethod fit
#
# @title "Fits the model used by this analyzer"
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
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fit", "GenomeAnalyzer", abstract=TRUE);





setConstructorS3("GenomeCnAnalyzer", function(..., mergeStrands=TRUE, combineAlleles=TRUE) {
  extend(GenomeAnalyzer(...), "GenomeCnAnalyzer",
    mergeStrands=mergeStrands, 
    combineAlleles=combineAlleles
  )
})



setConstructorS3("DChipCnAnalyzer", function(...) {
  extend(GenomeCnAnalyzer(...), "DChipCnAnalyzer",
    .plm = NULL
  )
})


setMethodS3("getModel", "DChipCnAnalyzer", function(this, ...) {
  # Check cache
  plm <- this$.plm;

  if (is.null(plm)) {
    ds <- getDataSet(this);
    plm <- MbeiCnPlm(ds, mergeStrands=this$mergeStrands, 
                         combineAlleles=this$combineAlleles);
    this$.plm <- plm;
  }

  plm;
})


setMethodS3("getChipEffects", "DChipCnAnalyzer", function(this, ...) {
  plm <- getModel(this);
  ces <- getChipEffects(plm, ...);
  ces;
})

setMethodS3("getProbeAffinities", "DChipCnAnalyzer", function(this, ...) {
  plm <- getModel(this);
  paf <- getProbeAffinities(plm, ...);
  paf;
})


setMethodS3("setup", "DChipCnAnalyzer", function(this, moreUnits=1, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Setting up ", class(this)[1], " object");

  model <- getModel(this);

  verbose && enter(verbose, "Get probe affinities");
  paf <- getProbeAffinities(this, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && enter(verbose, "Get chip effects");
  ces <- getChipEffects(this, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && enter(verbose, "Retrieving genome information");
  gi <- getGenomeInformation(this);
  verbose && exit(verbose);

  res <- list(model=model, paf=paf, ces=ces, gi=gi);

  verbose && exit(verbose);

  res;
}, protected=TRUE)




setMethodS3("fit", "DChipCnAnalyzer", function(this, chromosomes=getChromosomes(this), moreUnits=1, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosomes':
  chromosomes <- Arguments$getCharacters(chromosomes);
  chromosomes <- intersect(chromosomes, getChromosomes(this));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setup <- setup(this, verbose=verbose);
  attachLocally(setup);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit chromosome by chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(model);
  chipType <- getChipType(cdf);
  for (chr in chromosomes) {
    verbose && enter(verbose, "Chromosome ", chr);

    verbose && enter(verbose, "Identify units on chromosome");
    units <- getUnitIndices(gi, chromosome=chr);
    nunits <- length(units);
    verbose && printf(verbose, "Found %d units", nunits);
    verbose && str(verbose, units);
    verbose && exit(verbose);

    verbose && enter(verbose, sprintf("Fitting %d chromsome %s SNPs (%s)", 
                                                     nunits, chr, chipType));
    uDone <- fit(model, units=units, moreUnits=moreUnits, verbose=less(verbose));
    nDone <- length(uDone);
    verbose && printf(verbose, "Fitted %d SNPs\n", nDone);
    verbose && exit(verbose);

    if (nDone > 0) {
      verbose && enter(verbose, "Updating estimates of average chip effects");
      ceAvg <- getAverageFile(ces, units=uDone, force=TRUE, verbose=less(verbose));

      # Just in case for now. /HB 2006-10-31 (TODO)
      ceAvg$mergeStrands <- model$mergeStrands;
      ceAvg$combineAlleles <- model$combineAlleles;

      verbose && exit(verbose);
    }

    if (nDone > 0) {
      verbose && enter(verbose, "Updating \"units-todo\" cache");
      unitsTodo <- getUnitsTodo(this, verbose=less(verbose))[[chr]];
      unitsTodo$todo <- setdiff(unitsTodo$todo, uDone);
      this$.unitsTodo[[chr]] <- unitsTodo;
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Garbage collect");
    rm(units, nunits, uDone, nDone);
    gc();
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (chr in ...)
})


setMethodS3("fitGlad", "DChipCnAnalyzer", function(this, arrays=1:nbrOfArrays(this), chromosomes=getChromosomes(this), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosomes':
  chromosomes <- Arguments$getCharacters(chromosomes);
  chromosomes <- intersect(chromosomes, getChromosomes(this));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get chip effects
  ces <- getChipEffects(this, verbose=less(verbose));

  # Get reference
  ceAvg <- getAverageFile(ces, verbose=less(verbose));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Chromosome by chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list();
  for (chr in chromosomes) {
    res[[chr]] <- list();
    for (aa in seq(along=arrays)) {
      array <- arrays[aa];
      ce <- getFile(ces, array);
      verbose && enter(verbose, sprintf("Array %s (#%d of %d) on chromosome %s", 
                                      getName(ce), aa, length(arrays), chr));

      fit <- fitGlad(ce, reference=ceAvg, chromosome=chr, ..., verbose=less(verbose));

      verbose && enter(verbose, "Calling hooks");
      callHooks("fitGlad.DChipCnANalyzer.onFit", fit=fit, chromosome=chr, ce=ce);
      verbose && exit(verbose);
      
      res[[chr]][[aa]] <- fit;

      rm(fit, ce, array);

      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Garbage collect");
    gc();
    verbose && exit(verbose);
  } # for (chr in ...)
  res;
})


##############################################################################
# HISTORY:
# 2006-10-31
# o Created.
##############################################################################
