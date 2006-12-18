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
#   \item{tags}{A @character @vector of tags.}
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
setConstructorS3("GladModel", function(ces=NULL, reference=NULL, tags="*", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    tags[tags == "*"] <- "GLAD";
  }


  extend(Object(), "GladModel",
    .ces = ces,
    .reference = reference
  )
})

setMethodS3("as.character", "GladModel", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
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


setMethodS3("getName", "GladModel", function(this, ...) {
  # Get name of chip-effect set
  ces <- getChipEffects(this);
  getName(ces);
})


setMethodS3("getTags", "GladModel", function(this, ...) {
  ces <- getChipEffects(this);
  c(getTags(ces), this$.tags);
})


setMethodS3("getFullName", "GladModel", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "GladModel", function(this, ...) {
  "gladData";
})

setMethodS3("getPath", "GladModel", function(this, ...) {
  # Create the (sub-)directory tree for the dataset

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  cdf <- getCdf(this);
  chipType <- getChipType(cdf);
  chipType <- gsub("[,-]monocell$", "", chipType);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
})

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


setMethodS3("fit", "GladModel", function(this, arrays=1:nbrOfArrays(this), chromosomes=getChromosomes(this), ..., .retResults=TRUE, verbose=FALSE) {
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

  path <- getPath(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Chromosome by chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list();
  for (aa in seq(along=arrays)) {
    array <- arrays[aa];
    ce <- getFile(ces, array);
    arrayName <- getName(ce);
    res[[arrayName]] <- list();
    for (chr in chromosomes) {
      chrIdx <- match(chr, c(1:22, "X", "Y"));
      verbose && enter(verbose, 
                             sprintf("Array %s (#%d of %d) on chromosome %s", 
                                       arrayName, aa, length(arrays), chr));

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Get pathname
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Add tags chrNN,GLAD,<reference tags>
      tags <- setdiff(getTags(ce), "chipEffects");
      tags <- c(tags, sprintf("chr%02d", chrIdx));
      tags <- c(tags, "GLAD", getName(reference), getTags(reference));
      filename <- sprintf("%s.xdr", paste(c(arrayName, tags), collapse=","));
      pathname <- filePath(path, filename);

      # Already done?
      if (isFile(pathname)) {
        verbose && enter(verbose, "Loading results from file");
        verbose && cat(verbose, "Pathname: ", pathname);
        fit <- loadObject(pathname);
        verbose && exit(verbose);
      } else {
        fit <- fitGlad(ce, reference=reference, chromosome=chr, ...,
                                                      verbose=less(verbose));
        verbose && enter(verbose, "Saving to file");
        verbose && cat(verbose, "Pathname: ", pathname);
        saveObject(fit, file=pathname);
        verbose && exit(verbose);
      }

      verbose && enter(verbose, "Calling onFit() hooks");
      callHooks("onFit.fitGlad.GladModel", fit=fit, chromosome=chr, ce=ce, array=array);
      verbose && exit(verbose);

      if (.retResults)
        res[[arrayName]][[chr]] <- fit;

      rm(fit);

      verbose && exit(verbose);
    } # for (chr in ...)

    rm(ce, arrayName);
  } # for (aa in ...)

  invisible(res);
})


setMethodS3("plot", "GladModel", function(x, ..., pixelsPerMb=3, zooms=2^(0:7), pixelsPerTick=2.5, height=400, imageFormat="png", skip=TRUE, path=NULL, callSet=NULL, verbose=FALSE) {
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

  # Argument 'callSet':
  if (!is.null(callSet)) {
    if (!inherits(callSet, "GenotypeCallSet"))
      throw("Argument 'callSet' is not a GenotypeCallSet: ", class(callSet)[1]);
  }

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

  # Get chip type
  cdf <- getCdf(this);
  chipType <- getChipType(cdf);
  chipType <- gsub("[,-]monocell$", "", chipType);

  # The figure path
  if (is.null(path)) {
    path <- filePath(rootPath, getFullName(this), chipType, expandLinks="any");
    path <- filePath(path, expandLinks="any");
  }
  mkdirs(path);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the PNG device
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(imageFormat)) {
    imageFormat <- "current";
  }

  if (identical(imageFormat, "current")) {
    plotDev <- NULL;
  } else if (identical(imageFormat, "screen")) {
    screenDev <- function(pathname, width, height, ..., xpinch=50, ypinch=xpinch) {
      # Dimensions are in pixels. Rescale to inches
      width <- width/xpinch;
      height <- height/ypinch;
      x11(width=width, height=height, ...);
    }

    # When plotting to the screen, use only the first zoom
    zooms <- zooms[1];
    plotDev <- screenDev;
  } else if (identical(imageFormat, "png")) {
    pngDev <- System$findGraphicsDevice();
    plotDev <- pngDev;
  }


  callCols <- c("-"="lightgray", AA="red", AB="blue", BB="red", NC="orange");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Define the plot function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  on.exit({
    setHook("onFit.fitGlad.GladModel", NULL, action="replace");
  })

  setHook("onFit.fitGlad.GladModel", function(fit, chromosome, ce, array) {
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    tryCatch({

    # Get full name 
    fullname <- getFullName(ce);
    fullname <- gsub(",chipEffects$", "", fullname);

    chromosomeIdx <- match(chromosome, getChromosomes(this));
    nbrOfBases <- genome$nbrOfBases[chromosomeIdx];
    widthMb <- nbrOfBases / 1e6;

    verbose && enter(verbose, sprintf("Plotting %s for chromosome %s (%02d) [%.2fMb]", fullname, chromosome, chromosomeIdx, widthMb));

    for (zz in seq(along=zooms)) {
      zoom <- zooms[zz];

      # Create the pathname to the file
      imgName <- sprintf("%s,glad,chr%02d,x%04d.%s", fullname, chromosomeIdx, zoom, imageFormat);
      pathname <- filePath(path, imgName);

      # pngDev() (that is bitmap()) does not accept spaces in pathnames
      pathname <- gsub(" ", "_", pathname);
      if (!imageFormat %in% c("screen", "current")) {
        if (skip && isFile(pathname)) {
          next;
        }
      }
      # Calculate Mbs per ticks
      ticksBy <- 10^ceiling(log10(pixelsPerTick / (zoom * pixelsPerMb)));
      width <- as.integer(zoom * widthMb * pixelsPerMb);
      # Plot to PNG file
      verbose && printf(verbose, "Pathname: %s\n", pathname);
      verbose && printf(verbose, "Dimensions: %dx%d\n", width, height);
      verbose && printf(verbose, "Ticks by: %f\n", ticksBy);

      if (!is.null(plotDev))
        plotDev(pathname, width=width, height=height);

      tryCatch({
        verbose && enter(verbose, "Plotting graph");
        opar <- par(xaxs="r");
        plot(fit, ticksBy=ticksBy);
        stext(chipType, side=4, pos=1, line=0, cex=0.7, col="gray");
        verbose && enter(verbose, "Adding genotype calls");
        pv <- fit$profileValues;
        units <- pv$units;
        # Read the calls
        call <- callSet[units,array];  
        call <- as.character(call);

        x <- pv$PosBase; 
        x <- x/1e6;
        ylim <- par("usr")[3:4];
        ylim <- ylim + c(+1,-1)*0.04*diff(ylim);
        ylim <- ylim + c(+1,-1)*0.04*diff(ylim);
        y <- rep(ylim[1], length(callCols));
        names(y) <- names(callCols);
        y["AB"] <- y["AB"] + 0.02*diff(ylim);
        y <- y[call];
        points(x,y, pch=".", cex=2, col=callCols[call]);

        verbose && exit(verbose);
        verbose && exit(verbose);
      }, error = function(ex) {
        print(ex);
      }, finally = {
        par(opar);
        if (!imageFormat %in% c("screen", "current"))
          dev.off();
      });
    } # for (zoom ...)

    verbose && exit(verbose);
    }, error = function(ex) {
      cat("ERROR caught in onFit.fitGlad.GladModel():\n");
      print(ex);
    }) # tryCatch()
  }, action="replace")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Start fitting and plotting
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit(this, ..., .retResults=FALSE, verbose=verbose);

  invisible();
})


setMethodS3("getLog2Ratios", "GladModel", function(this, ..., verbose=FALSE) {
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
  })

  res;
}) # getLog2Ratios()



setMethodS3("getRegions", "GladModel", function(this, ..., url="ucsc", margin=10, flat=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'url':
  if (identical(url, "ucsc")) {
    url <- function(chromosome, start, stop) {
      uri <- "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=vertebrate&org=Human&db=hg18";
      sprintf("%s&position=chr%s%%3A%d-%d", uri, chromosome, as.integer(start), as.integer(stop));
    }
  }

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

    # Add URL?
    if (!is.null(url)) {
      chromosome <- df[,"Chromosome"];
      start <- df[,"start"];
      stop <- df[,"stop"];
      m <- margin*abs(stop-start);
      start <- start-m;
      start[start < 0] <- 0;
      stop <- stop + m;
      urls <- character(nrow(df));
      for (rr in seq(along=urls)) { 
        urls[rr] <- url(chromosome[rr], start[rr], stop[rr]);
      }
      df <- cbind(df, url=urls);
    }

    df;
  })

  if (flat) {
    df <- NULL;
    for (kk in seq(along=res)) {
      df <- rbind(df, cbind(sample=names(res)[kk], res[[kk]]));
      res[[kk]] <- NA;
    }
    row.names(df) <- seq(length=nrow(df));
    res <- df;
  }  

  res;
})


setMethodS3("writeRegions", "GladModel", function(this, arrays=1:nbrOfArrays(this), format=c("xls", "wig"), digits=3, ..., oneFile=FALSE, skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  arrays <- Arguments$getIndices(arrays, range=c(1,nbrOfArrays(this)));

  # Argument 'format':
  format <- match.arg(format);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Add a "GLAD" tag.
  dataSetName <- getFullName(this);

  ces <- getChipEffects(this);
  arrayNames <- getNames(ces);

  path <- getPath(this);
  mkdirs(path);

  if (oneFile) {
    filename <- sprintf("%s,GLAD,regions.%s", getFullName(this), format); 
    pathname <- filePath(path, filename);
    pathname <- Arguments$getWritablePathname(pathname);
    if (!skip && isFile(pathname))
      file.remove(pathname);
  }

  res <- list();
  for (array in arrays) {
    name <- arrayNames[array];
    verbose && enter(verbose, sprintf("Array #%d (of %d) - %s", 
                                             array, length(arrays), name));
    df <- getRegions(this, arrays=array, ..., verbose=less(verbose))[[1]];
    names(df) <- gsub("Smoothing", "log2CN", names(df));

    if (identical(format, "xls")) {
      # Append column of sample names?
      if (oneFile)
        df <- cbind(sample=name, df);
    } else if (identical(format, "wig")) {
      # Write a four column WIG/BED table
      df <- df[,c("Chromosome", "start", "stop", "log2CN")];

      # In the UCSC Genome Browser, the maximum length of one element
      # is 10,000,000 bases.  Chop up long regions in shorter contigs.
      verbose && enter(verbose, sprintf("Chopping up too long segment"));
      MAX.LENGTH = 10e6-1;
      start <- df[,"start"];
      stop <- df[,"stop"];
      len <- stop-start;
      tooLong <- which(len > MAX.LENGTH);
      if (length(tooLong) > 0) {
        dfXtra <- NULL;
        for (rr in tooLong) {
          x0 <- start[rr];
          while (x0 < stop[rr]) {
            x1 <- min(x0 + MAX.LENGTH, stop[rr]);
            df1 <- df[rr,];
            df1[,"start"] <- x0;
            df1[,"stop"] <- x1;
            dfXtra <- rbind(dfXtra, df1);          
            x0 <- x1+1;
          }
        }
        df <- df[-tooLong,];
        df <- rbind(df, dfXtra);
        rm(dfXtra);
        row.names(df) <- seq(length=nrow(df));
      }
      verbose && exit(verbose);
      # Make sure the items are ordered correctly
      chrIdx <- as.integer(df[,"Chromosome"]);
      o <- order(chrIdx, df[,"start"]);
      df <- df[o,];

      # All chromosomes should have prefix 'chr'.
      chrIdx <- as.integer(df[,"Chromosome"]);
      df[chrIdx == 23,"Chromosome"] <- "X";
      df[,"Chromosome"] <- paste("chr", df[,"Chromosome"], sep="");
    }

    # Apply digits
    for (cc in seq(length=ncol(df))) {
      value <- df[,cc];
      if (is.double(value)) {
        df[,cc] <- round(value, digits=digits);
      }
    }

    if (!oneFile) {
      savename <- name;
      filename <- sprintf("%s,GLAD,regions.%s", savename, format); 
      pathname <- filePath(path, filename);
      if (!oneFile && !skip && isFile(pathname))
        file.remove(pathname);
    }

    # Writing to file
    verbose && cat(verbose, "Pathname: ", pathname);
    if (identical(format, "xls")) {
      col.names <- (array == arrays[1]);
      write.table(df, file=pathname, sep="\t", col.names=col.names, row.names=FALSE, quote=FALSE, append=oneFile);
    } else if (identical(format, "wig")) {
      # Write track control
      trackAttr <- c(type="wiggle_0");
      trackAttr <- c(trackAttr, name=sprintf("\"%s\"", name));
      trackAttr <- c(trackAttr, group="\"GLAD regions\"");
      trackAttr <- c(trackAttr, priority=array);
      trackAttr <- c(trackAttr, graphType="bar");
      trackAttr <- c(trackAttr, visibility="full");
      trackAttr <- c(trackAttr, maxHeightPixels="128:96:64");
      trackAttr <- c(trackAttr, yLineOnOff="on");
# HARD WIRED FOR NOW.  TO DO /hb 2006-11-27
col <- c("117,112,179", "231,41,138");
ylim <- c(-1,1);
      if (!is.null(col)) {
        trackAttr <- c(trackAttr, color=col[1], altColor=col[2]);
      }
      if (!is.null(ylim)) {
        trackAttr <- c(trackAttr, autoScale="off", 
              viewLimits=sprintf("%.2f:%.2f ", ylim[1], ylim[2]));
      }
      trackAttr <- paste(names(trackAttr), trackAttr, sep="=");
      trackAttr <- paste(trackAttr, collapse=" ");
      trackAttr <- paste("track ", trackAttr, "\n", sep="");
      verbose && str(verbose, trackAttr);
      cat(file=pathname, trackAttr, append=oneFile);

      # Write data
      verbose && str(verbose, df);
      write.table(df, file=pathname, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE, append=oneFile);
    }
    verbose && exit(verbose);
    res[[array]] <- df;
  }

  invisible(res);
})


setMethodS3("writeWig", "GladModel", function(this, ...) {
  ces <- getChipEffects(this);
  reference <- getReference(this);
  writeWig(ces, reference=reference, ...);
}) # writeWig()



##############################################################################
# HISTORY:
# 2006-11-29
# o Added chip type annotation to plot() and option to plot to screen.
# 2006-11-27
# o Added argument 'flat' to getRegions().
# 2006-11-24
# o Now the fit() function of GladModel stores the profileCGH object as a
#   binary XDR file in the default path, see getPath().
# 2006-11-23
# o Added writeWig().
# 2006-11-22
# o Added writeRegions().
# o Added fit(), plot(), and getRegions().
# o Re-created from the CnAnalyzer class from 2006-10-31.
##############################################################################
