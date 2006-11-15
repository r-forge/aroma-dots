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
    .dataSets = list(raw=dataSet)
  )
})

setMethodS3("as.character", "GenomeAnalyzer", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  s <- c(s, "DATASET:");
  s <- c(s, as.character(getDataSet(this)));
  class(s) <- "GenericSummary";
  s;
})

setMethodS3("getDataSets", "GenomeAnalyzer", function(this, ...) {
  this$.dataSets;
}, protected=TRUE)

setMethodS3("appendDataSet", "GenomeAnalyzer", function(this, ...) {
  this$.dataSets <- c(this$.dataSets, list(...));
}, protected=TRUE)

setMethodS3("getDataSet", "GenomeAnalyzer", function(this, ...) {
  # Gets the last data set
  dss <- getDataSets(this);
  dss[[length(dss)]];
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
setMethodS3("getUnitsTodo", "GenomeAnalyzer", function(this, chromosomes=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'chromosomes':
  if (is.null(chromosomes)) {
    chromosomes <- getChromosomes(this);
  } else {
    missing <- !(chromosomes %in% getChromosomes(this));
    if (any(missing)) {
      throw("Argument 'chromosomes' contains unknown elements: ", 
                                 paste(chromosomes[missing], collapse=", "));
    }
  }

  unitsTodo <- this$.unitsTodo;
  if (force || is.null(unitsTodo) || !is.null(chromosomes)) {
    model <- getModel(this);
    unitsTodo <- list();
    gi <- getGenomeInformation(this);
    for (chr in chromosomes) {
      verbose && enter(verbose, "Chromosome ", chr);
  
      units <- getUnitIndices(gi, chromosome=chr);
      verbose && printf(verbose, "Found in total %d units on chromosome", length(units));
      verbose && str(verbose, units);
  
      utodo <- findUnitsTodo(model, units=units, verbose=less(verbose));
      verbose && printf(verbose, "Out of which %d units are not fitted.\n", length(utodo));
  
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




setConstructorS3("CnPlmAnalyzer", function(...) {
  extend(GenomeCnAnalyzer(...), "CnPlmAnalyzer",
    .plm = NULL
  )
})

setMethodS3("getChipEffects", "CnPlmAnalyzer", function(this, ...) {
  plm <- getModel(this);
  ces <- getChipEffects(plm, ...);
  ces;
})

setMethodS3("getProbeAffinities", "CnPlmAnalyzer", function(this, ...) {
  plm <- getModel(this);
  paf <- getProbeAffinities(plm, ...);
  paf;
})


setMethodS3("setup", "CnPlmAnalyzer", function(this, moreUnits=1, ..., verbose=FALSE) {
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



setMethodS3("fit", "CnPlmAnalyzer", function(this, chromosomes=getChromosomes(this), moreUnits=1, ..., verbose=FALSE) {
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
  # Identify units todo (cached)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying units to be estimated");
  getUnitsTodo(this, verbose=less(verbose));
  verbose && exit(verbose);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit chromosome by chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(model);
  chipType <- getChipType(cdf);
  for (chr in chromosomes) {
    verbose && enter(verbose, "Chromosome ", chr);

    unitsInfo <- getUnitsTodo(this, verbose=less(verbose))[[chr]];
    nAllUnits <- length(unitsInfo$units);
    units <- unitsInfo$todo;
    nunits <- length(units);
    verbose && printf(verbose, "Found %d out of %d units on chromosome %s to fit.\n", nunits, nAllUnits, chr);
    if (nunits > 0)
      verbose && str(verbose, units);

    if (nunits > 0) {
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
      rm(uDone, nDone);
      gc();
      verbose && exit(verbose);
    } # if (nunits > 0)
    rm(units, nunits);

    verbose && exit(verbose);
  } # for (chr in ...)
})


setMethodS3("fitGlad", "CnPlmAnalyzer", function(this, arrays=1:nbrOfArrays(this), chromosomes=getChromosomes(this), ..., verbose=FALSE) {
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
      callHooks("fitGlad.CnPlmAnalyzer.onFit", fit=fit, chromosome=chr, ce=ce);
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


setMethodS3("plotGlad", "CnPlmAnalyzer", function(this, pixelsPerMb=3, zooms=2^(0:8), pixelsPerTick=2.5, height=400, ..., skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pixelsPerMb':
  pixelsPerMb <- Arguments$getDouble(pixelsPerMb, range=c(0.001,9999));

  # Argument 'zooms':
  zooms <- Arguments$getIntegers(zooms, range=c(1,9999));
  zooms <- unique(zooms);

  # Argument 'pixelsPerMb':
  pixelsPerTick <- Arguments$getDouble(pixelsPerTick, range=c(1,256));

  # Argument 'height':
  height <- Arguments$getInteger(height, range=c(1,4096));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get genome annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  genome <- readTable("annotations/hgChromosomes.txt", header=TRUE, 
                            colClasses=c(nbrOfBases="integer"), row.names=1);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Output path
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rootPath <- "chromosomeExplorer";

  # Get name of chip-effect set
  ces <- getChipEffects(this, verbose=less(verbose));
  fullname <- getFullName(ces);

  # Get chip type
  cdf <- getCdf(ces);
  chipType <- getChipType(cdf);

  # The glad directory
  chipType <- gsub("-monocell", "", chipType);
  path <- filePath(rootPath, fullname, chipType);
  mkdirs(path);

  # The figure path
  figPath <- filePath(path);
  mkdirs(figPath);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the PNG device
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pngDev <- System$findGraphicsDevice();
  imgExt <- "png";


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Define the plot function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  on.exit({
    setHook("fitGlad.CnPlmAnalyzer.onFit", NULL, action="replace");
  })

  setHook("fitGlad.CnPlmAnalyzer.onFit", function(fit, chromosome, ce) {
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    # Get full name 
    fullname <- getFullName(ce);
    fullname <- gsub(",chipEffects$", "", fullname);

    chromosomeIdx <- match(chromosome, getChromosomes(this));
    nbrOfBases <- genome$nbrOfBases[chromosomeIdx];
    widthMb <- nbrOfBases / 1e6;

    verbose && enter(verbose, sprintf("Plotting %s for chromosome %s (%02d) [%.2fMb]", fullname, chromosome, chromosomeIdx, widthMb));

    for (zz in seq(along=zooms)) {
      zoom <- zooms[zz];
      imgName <- sprintf("%s,glad,chr%02d,x%04d.%s", fullname, chromosomeIdx, zoom, imgExt);
      pathname <- filePath(figPath, imgName);

      # pngDev() (that is bitmap()) does not accept spaces in pathnames
      pathname <- gsub(" ", "_", pathname);

      if (skip && isFile(pathname)) {
        next;
      }

      # Calculate Mbs per ticks
      ticksBy <- 10^ceiling(log10(pixelsPerTick / (zoom * pixelsPerMb)));
      width <- as.integer(zoom * widthMb * pixelsPerMb);
      # Plot to PNG file
      verbose && printf(verbose, "Pathname: %s\n", pathname);
      verbose && printf(verbose, "Dimensions: %dx%d\n", width, height);
      verbose && printf(verbose, "Ticks by: %f\n", ticksBy);

      pngDev(pathname, width=width, height=height);
      tryCatch({
        verbose && enter(verbose, "Plotting graph");
        opar <- par(xaxs="r");
        plot(fit, ticksBy=ticksBy);
        verbose && exit(verbose);
      }, error = function(ex) {
         print(ex);
      }, finally = {
         par(opar);
        dev.off();
      });
    } # for (zoom ...)
    verbose && exit(verbose);
  }, action="replace")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Start fitting and plotting
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- fitGlad(this, ..., verbose=verbose);

  invisible(res);
})




setConstructorS3("MbeiCnAnalyzer", function(...) {
  extend(CnPlmAnalyzer(...), "MbeiCnAnalyzer")
})

setMethodS3("getModel", "MbeiCnAnalyzer", function(this, ...) {
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

setMethodS3("preprocess", "MbeiCnAnalyzer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  ds <- getDataSet(this);
  norm <- QuantileNormalizer(ds, ..., subsetToAvg=1/3);
  dsN <- process(norm, verbose=less(verbose));
  dsN$normalizer <- norm;
  print(dsN);

  # Append this data set
  appendDataSet(this, preprocess=dsN);

  invisible(dsN);
})



setConstructorS3("RmaCnAnalyzer", function(...) {
  extend(CnPlmAnalyzer(...), "RmaCnAnalyzer")
})

setMethodS3("getModel", "RmaCnAnalyzer", function(this, ...) {
  # Check cache
  plm <- this$.plm;

  if (is.null(plm)) {
    ds <- getDataSet(this);
    plm <- RmaCnPlm(ds, mergeStrands=this$mergeStrands, 
                         combineAlleles=this$combineAlleles);
    this$.plm <- plm;
  }

  plm;
})


setMethodS3("preprocess", "RmaCnAnalyzer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  ds <- getDataSet(this);
  norm <- QuantileNormalizer(ds, ..., subsetToAvg=1/3);
  dsN <- process(norm, verbose=less(verbose));
  dsN$normalizer <- norm;
  print(dsN);

  # Append this data set
  appendDataSet(this, preprocess=dsN);

  invisible(dsN);
})




##############################################################################
# HISTORY:
# 2006-11-07
# o Now plotGlad() removed the hook function when exiting.
# o Now MbeiCnAnalyzer and RmaCnAnalyzer subclasses CnPlmAnalyzer.
# 2006-10-31
# o Created.
##############################################################################
