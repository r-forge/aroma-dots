###########################################################################/**
# @RdocClass GladModel
#
# @title "The GladModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the GLAD [1] model.
# }
# 
# @synopsis
#
# \arguments{
#   \item{ces}{A @see "CnChipEffectSet".}
#   \item{reference}{A @see "CnChipEffectFile".}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \references{
#  [1] Hupe P et al. \emph{Analysis of array CGH data: from signal ratio to
#      gain and loss of DNA regions}. Bioinformatics, 2004, 20, 3413-3422.\cr
# }
#
#*/###########################################################################
setConstructorS3("GladModel", function(ces=NULL, reference=NULL, ...) {
  # Argument 'ces':
  if (!is.null(ces)) {
    if (!inherits(ces, "CnChipEffectSet")) {
      throw("Argument 'ces' is not a CnChipEffectSet object: ", class(ces));
    }

    # Currently only total copy-number estimates are accepted
    if (!ces$combineAlleles) {
      throw("Unsupported copy-number chip effects. Currently only total copy-number estimates are supported: ces$combineAlleles == FALSE");
    }
  }

  # Argument 'reference':
  if (!is.null(reference)) {
    if (!inherits(reference, "CnChipEffects")) {
      throw("Argument 'reference' is not a CnChipEffects object: ", 
                                                            class(reference));
    }

    if (reference$combineAlleles != ces$combineAlleles) {
       throw("The reference chip effects are not compatible with the chip-effect set. One is combining the alleles the other is not.");
    }

    if (reference$mergeStrands != ces$mergeStrands) {
       throw("The reference chip effects are not compatible with the chip-effect set. One is merging the strands the other is not.");
    }
  }


  extend(Object(), "GladModel",
    .ces = ces,
    .reference = reference
  )
})

setMethodS3("as.character", "GladModel", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, "Chip effects:");
  s <- c(s, as.character(getChipEffects(this)));
  s <- c(s, "Reference:");
  if (is.null(getReference(this))) {
    s <- c(s, "<average across arrays>");
  } else {
    s <- c(s, as.character(getReference(this)));
  }
  s <- c(s, "Genome information:", as.character(getGenomeInformation(this)));
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
})

setMethodS3("getChipEffects", "GladModel", function(this, ...) {
  this$.ces;
})

setMethodS3("getReference", "GladModel", function(this, ...) {
  this$.reference;
})

setMethodS3("setReference", "GladModel", function(this, reference=NULL, ...) {
  # Argument 'reference':
  if (!is.null(reference)) {
    if (!inherits(reference, "CnChipEffects")) {
      throw("Argument 'reference' is not a CnChipEffects object: ", 
                                                            class(reference));
    }
  }

  this$.reference <- reference;
})

setMethodS3("nbrOfArrays", "GladModel", function(this, ...) {
  ces <- getChipEffects(this);
  nbrOfArrays(ces);
})

setMethodS3("getCdf", "GladModel", function(this, ...) {
  ces <- getChipEffects(this);
  getCdf(ces);
}, protected=TRUE)

setMethodS3("getChipType", "GladModel", function(this, ...) {
  cdf <- getCdf(this);
  getChipType(cdf);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getChromosomes
#
# @title "Gets the chromosomes available"
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
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChromosomes", "GladModel", function(static, ...) {
  c(1:22, "X");
}, static=TRUE)



setMethodS3("getGenomeInformation", "GladModel", function(this, ...) {
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
  cdf <- getCdf(this);
  gi <- getGenomeInformation(cdf);
  verbose && exit(verbose);

  gi;
})


setMethodS3("fit", "GladModel", function(this, arrays=1:nbrOfArrays(this), chromosomes=getChromosomes(this), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  arrays <- Arguments$getIndices(arrays, range=c(1,nbrOfArrays(this)));

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
  # Retrieving chip effects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get chip effects
  ces <- getChipEffects(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving reference chip effects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  reference <- getReference(this);
  if (is.null(reference)) {
    verbose && enter(verbose, "No reference specified. Calculating average chip effects");
    reference <- getAverageFile(ces, verbose=less(verbose));
    verbose && exit(verbose);
  }
  verbose && cat(verbose, "Using reference:");
  verbose && print(verbose, reference);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Chromosome by chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list();
  for (aa in seq(along=arrays)) {
    ce <- getFile(ces, arrays[aa]);
    arrayName <- getName(ce);
    res[[arrayName]] <- list();
    for (chr in chromosomes) {
      verbose && enter(verbose, 
                             sprintf("Array %s (#%d of %d) on chromosome %s", 
                                       arrayName, aa, length(arrays), chr));

      fit <- fitGlad(ce, reference=reference, chromosome=chr, ...,
                                                      verbose=less(verbose));
      res[[arrayName]][[chr]] <- fit;

      verbose && enter(verbose, "Calling onFit() hooks");
      callHooks("onFit.fitGlad.GladModel", fit=fit, chromosome=chr, ce=ce);
      verbose && exit(verbose);

      rm(fit);

      verbose && exit(verbose);
    } # for (chr in ...)

    rm(ce, arrayName);
  } # for (aa in ...)

  invisible(res);
})


setMethodS3("plot", "GladModel", function(x, ..., pixelsPerMb=3, zooms=2^(0:7), pixelsPerTick=2.5, height=400, skip=TRUE, verbose=FALSE) {
  # To please R CMD check.
  this <- x;

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
  ces <- getChipEffects(this);
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
    setHook("onFit.fitGlad.GladModel", NULL, action="replace");
  })

  setHook("onFit.fitGlad.GladModel", function(fit, chromosome, ce) {
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
  res <- fit(this, ..., verbose=verbose);

  invisible(res);
})


setMethodS3("getRegions", "GladModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the regions for each of the GLAD fits (per array)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  suppressWarnings({
    res <- fit(this, ..., verbose=less(verbose,10));
  })
  res <- lapply(res, FUN=function(arrayFits) {
    df <- NULL;
    for (fit in arrayFits) {
      suppressWarnings({
        df0 <- getRegions(fit, ...);
      })
      df <- rbind(df, df0);
    }
    rownames(df) <- seq(length=nrow(df));
    df;
  })

  res;
})


setMethodS3("writeRegions", "GladModel", function(this, arrays=1:nbrOfArrays(this), ext="xls", digits=3, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  arrays <- Arguments$getIndices(arrays, range=c(1,nbrOfArrays(this)));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  ces <- getChipEffects(this);
  dataSetName <- getFullName(ces);
  # Add a "GLAD" tag.
  dataSetName <- sprintf("%s,GLAD", dataSetName);
  cdf <- getCdf(ces);
  chipType <- getChipType(cdf);
  chipType <- gsub("-monocell", "", chipType);
  arrayNames <- getNames(ces);

  path <- filePath("glad", dataSetName, chipType, expandLinks="any");
  mkdirs(path);

  res <- list();
  for (array in arrays) {
    name <- arrayNames[array];
    verbose && enter(verbose, sprintf("Array #%d (of %d) - %s", 
                                             array, length(arrays), name));
    df <- getRegions(this, arrays=array, ..., verbose=less(verbose))[[1]];
    names(df) <- gsub("Smoothing", "log2CN", names(df));

    # Apply digits
    for (cc in seq(length=ncol(df))) {
      value <- df[,cc];
      if (is.double(value)) {
        df[,cc] <- round(value, digits=digits);
      }
    }

    # Writing to file
    filename <- sprintf("%s,GLAD,regions.%s", name, ext); 
    pathname <- filePath(path, filename);
    verbose && cat(verbose, "Pathname: ", pathname);
    write.table(df, file=pathname, sep="\t", row.names=FALSE, quote=FALSE);
    verbose && exit(verbose);
    res[[array]] <- df;
  }

  invisible(res);
})



##############################################################################
# HISTORY:
# 2006-11-22
# o Added writeRegions().
# o Added fit(), plot(), and getRegions().
# o Re-created from the CnAnalyzer class from 2006-10-31.
##############################################################################
