###########################################################################/**
# @RdocClass GladModel
#
# @title "The GladModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Gain and Loss Analysis of DNA regions 
#  (GLAD) model [1].
#  This class can model chip-effect estimates obtained from multiple
#  chip types, and not all samples have to be available on all chip types.
# }
# 
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "ChipEffectSetTuple".}
#   \item{...}{Arguments passed to the constructor of 
#              @see "CopyNumberSegmentationModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   Data from multiple chip types are combined "as is".  This is based
#   on the assumption that the relative chip effect estimates are
#   non-biased (or at the equally biased across chip types).
#   Note that in GLAD there is no way to down weight certain data points,
#   which is why we can control for differences in variance across
#   chip types.
# }
#
# \section{Benchmarking}{
#   In high-density copy numbers analysis, the most time consuming step
#   is fitting the GLAD model.  The complexity of the model grows 
#   more than linearly (squared? exponentially?) with the number of data
#   points in the chromosome and sample being fitted.  This is why it
#   take much more than twice the time to fit two chip types together 
#   than separately.
# }
#
# @author
# 
# \seealso{
#  @see "CopyNumberSegmentationModel".
# }
#
# \references{
#  [1] Hupe P et al. \emph{Analysis of array CGH data: from signal ratio to
#      gain and loss of DNA regions}. Bioinformatics, 2004, 20, 3413-3422.\cr
# }
#*/###########################################################################
setConstructorS3("GladModel", function(cesTuple=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(cesTuple)) {
    require("GLAD") || throw("Package not loaded: GLAD");
  }

  extend(CopyNumberSegmentationModel(cesTuple=cesTuple, ...), "GladModel")
})


setMethodS3("getAsteriskTag", "GladModel", function(this, ...) {
  "GLAD";
}, protected=TRUE)


setMethodS3("plot", "GladModel", function(x, ..., pixelsPerMb=3, zooms=2^(0:7), pixelsPerTick=2.5, height=400, xmargin=c(50,50), imageFormat="current", skip=TRUE, path=NULL, callList=NULL, verbose=FALSE) {
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

  # Argument 'callList':
  chipTypes <- getChipTypes(this);
  if (length(callList) > 0) {
    if (!is.list(callList))
      callList <- list(callList);

    if (length(callList) != length(chipTypes)) {
      throw("Number of elements in argument 'callList' does not match the number of chip types: ", length(callList), " != ", length(chipTypes));
    }

    if (is.null(names(callList)))
      names(callList) <- chipTypes;

    for (chipType in chipTypes) {
      callSet <- callList[[chipType]];
      if (!is.null(callSet)) {
        if (!inherits(callSet, "GenotypeCallSet"))
          throw("Argument 'callList' contains a non-GenotypeCallSet: ", class(callSet)[1]);

        if (getChipType(callSet) != chipType) {
          throw("Argument 'callList' contains a GenotypeCallSet for a different chip type than the corresponding chip-effect set: ", getChipType(callSet), " != ", chipType);
        }
      }
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Get genome annotation data (chromosome lengths etc)
  genome <- getGenomeData(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Output path
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The report path
  if (is.null(path))
    path <- getReportPath(this);
  mkdirs(path);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the PNG device
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(imageFormat)) {
    imageFormat <- "current";
  }

  resScale <- 1;
  if (identical(imageFormat, "current")) {
    plotDev <- NULL;
    zooms <- zooms[1];
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
    if (identical(pngDev, png2))
      resScale <- 2;
  }


  callCols <- c("-"="lightgray", AA="red", AB="blue", BB="red", NC="orange");

  # Get chip type (used to annotate the plot)
  chipType <- getChipType(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Define the plot function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  on.exit({
    setHook("onFit.GladModel", NULL, action="replace");
  })

  setHook("onFit.GladModel", function(fit, fullname) {
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    tryCatch({
      # Extract the array name from the full name
      arrayFullName <- gsub("^(.*),chr[0-9][0-9].*$", "\\1", fullname);
      arrayName <- gsub("^([^,]*).*$", "\\1", arrayFullName);
  
      # Extract the chromosome from the GLAD fit object
      pv <- fit$profileValues;
      chromosome <- unique(pv$Chromosome);  # Should only be one!
##      chromosome[chromosome == "23"] <- "X";  # TODO
  
      # Infer the length (in bases) of the chromosome
##      chromosomeIdx <- match(chromosome, getChromosomes(this));
      nbrOfBases <- genome$nbrOfBases[chromosome];
      widthMb <- nbrOfBases / 1e6;
  
      verbose && enter(verbose, sprintf("Plotting %s for chromosome %02d [%.2fMB]", arrayName, chromosome, widthMb));
  
      for (zz in seq(along=zooms)) {
        zoom <- zooms[zz];
  
        # Create the pathname to the file
        imgName <- sprintf("%s,chr%02d,x%04d.%s", 
                          arrayFullName, chromosome, zoom, imageFormat);
        pathname <- filePath(path, imgName);
  
        # pngDev() (that is bitmap()) does not accept spaces in pathnames
        pathname <- gsub(" ", "_", pathname);
        if (!imageFormat %in% c("screen", "current")) {
          if (skip && isFile(pathname)) {
            next;
          }
        }
        # Calculate MBs per ticks
        ticksBy <- 10^ceiling(log10(pixelsPerTick / (zoom * pixelsPerMb)));
        xlimDiff <- (zoom * widthMb * pixelsPerMb);
        width <- round(zoom * widthMb * pixelsPerMb + sum(xmargin));
        # Plot to PNG file
        verbose && printf(verbose, "Pathname: %s\n", pathname);
        verbose && printf(verbose, "Dimensions: %dx%d\n", width, height);
        verbose && printf(verbose, "Ticks by: %f\n", ticksBy);
  
        if (!is.null(plotDev))
          plotDev(pathname, width=width, height=height);
        tryCatch({
          verbose && enter(verbose, "Plotting graph");
          opar <- par(xaxs="r");
          suppressWarnings({
            plot(fit, ticksBy=ticksBy, ..., xmargin=xmargin, resScale=resScale, flavor="ce");
          });
          stext(chipType, side=4, pos=1, line=0, cex=0.8);
          if (!is.null(callList)) {
            verbose && enter(verbose, "Adding genotype calls");
  
            # Get (chip type, unit) information
            chipType <- pv$chipType;
            unit <- pv$unit;
  
            # Figure out where to put the genotype track
            ylim <- par("usr")[3:4];
            ylim <- ylim + c(+1,-1)*0.04*diff(ylim);
            ylim <- ylim + c(+1,-1)*0.04*diff(ylim);
  
            for (chipType in chipTypes) {
              # Got genotype calls for this chip type?
              callSet <- callList[[chipType]];
              if (is.null(callSet))
                next;
  
              # Got chip-effect estimates for this chip type?
              idxs <- which(pv$chipType == chipType);
              if (length(idxs) == 0)
                next;
  
              # Got genotype calls for this array?
              if (!arrayName %in% getNames(callSet))
                next;
  
              # Get subset of genotype calls for this array & chromosome.
              units <- pv$unit[idxs];
              call <- callSet[units, arrayName];
              call <- as.character(call);
              # Get the positions of these calls
              x <- pv$PosBase[idxs]; 
              x <- x/1e6;
  
              # Plot the homozygote/heterozygote genotype tracks
              y <- rep(ylim[1], length(callCols));
              names(y) <- names(callCols);
              y["AB"] <- y["AB"] + 0.02*diff(ylim);
              y <- y[call];
              points(x,y, pch=".", cex=2, col=callCols[call]);
   
              rm(idxs, call, callSet, units, x, y);  # Not needed anymore
            } # for (chipType ...)
            verbose && exit(verbose);
          } # if (!is.null(callList))
  
          verbose && exit(verbose);
        }, error = function(ex) {
          print(ex);
        }, finally = {
          par(opar);
          if (!imageFormat %in% c("screen", "current"))
            dev.off();
        });
      } # for (zz in ...)
    }, error = function(ex) {
      cat("ERROR caught in onFit.GladModel():\n");
      print(ex);
    }, finally = {
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);
    }) # tryCatch()
  }, action="replace")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Start fitting and plotting
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit(this, ..., .retResults=FALSE, verbose=verbose);

  invisible();
})


setMethodS3("getRegions", "GladModel", function(this, ..., url="ucsc", organism="Human", hgVersion="hg17", margin=10, flat=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'url':
  if (identical(url, "ucsc")) {
    # The UCSC browser accepts chromsomes either 'X' or 23.  In other words,
    # we can stick with integers to be more general.
    url <- function(chromosome, start, stop) {
      uri <- "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=vertebrate";
      sprintf("%s&org=%s&db=%s&position=chr%s%%3A%d-%d", uri, organism, hgVersion, chromosome, as.integer(start), as.integer(stop));
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
  verbose && enter(verbose, "Obtaining GLAD fits (or fit if missing)");
  suppressWarnings({
    res <- fit(this, ..., .retResults=TRUE, verbose=less(verbose,10));
  })
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting regions from all fits");
  res <- lapply(res, FUN=function(arrayFits) {
    df <- NULL;
    # For each chromosome
    for (kk in seq(along=arrayFits)) {
      fit <- arrayFits[[kk]];
      if (!is.null(fit)) {
        verbose && enter(verbose, "Extracting regions for chromosome #", kk);
        suppressWarnings({
          df0 <- getRegions(fit, ...);
        })
        df <- rbind(df, df0);
        verbose && exit(verbose);
      }
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
  verbose && exit(verbose);

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


setMethodS3("writeRegions", "GladModel", function(this, arrays=NULL, format=c("xls", "wig"), digits=3, ..., oneFile=TRUE, skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  arrays <- indexOfArrays(this, arrays=arrays);

  # Argument 'format':
  format <- match.arg(format);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Setup
  fullname <- getFullName(this);
  arrayNames <- getNames(this);

  path <- getPath(this);
  mkdirs(path);

  if (oneFile) {
    filename <- sprintf("%s,regions.%s", fullname, format); 
    pathname <- filePath(path, filename);
    pathname <- Arguments$getWritablePathname(pathname);
    if (!skip && isFile(pathname)) {
      file.remove(pathname);
    }
  }

  res <- list();
  for (aa in seq(along=arrays)) {
    array <- arrays[aa];
    name <- arrayNames[array];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                               aa, name, length(arrays)));
    df <- getRegions(this, arrays=array, ..., verbose=less(verbose))[[1]];
    names(df) <- gsub("Smoothing", "log2", names(df));

    if (nrow(df) > 0) {
      if (identical(format, "xls")) {
        # Append column with sample names
        df <- cbind(sample=name, df);
      } else if (identical(format, "wig")) {
        # Write a four column WIG/BED table
        df <- df[,c("Chromosome", "start", "stop", "log2")];
  
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
        ## df[chrIdx == 23,"Chromosome"] <- "X"; ## REMOVED 2007-03-15
        df[,"Chromosome"] <- paste("chr", df[,"Chromosome"], sep="");
      }
  
      # Apply digits
      for (cc in seq(length=ncol(df))) {
        value <- df[,cc];
        if (is.double(value)) {
          df[,cc] <- round(value, digits=digits);
        }
      }
    } # if (nrow(df) > 0)

    if (!oneFile) {
      savename <- name;
      filename <- sprintf("%s,regions.%s", savename, format); 
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

  invisible(pathname);
})



##############################################################################
# HISTORY:
# 2007-08-20
# o Initial tests show that the updated GladModel gives identical results.
# o Now GladModel inherits from CopyNumberSegmentationModel.
# 2007-05-10
# o BUG FIX: getRegions() and getLog2Ratios() would give an error if a subset
#   of the chromosomes where queried.
# o Added more verbose output to getRegions().
# 2007-04-12
# o Now plot() of the GladModel writes the chip type annotation in black and
#   not in gray as before.
# 2007-03-24
# o BUG FIX: getPath() created the root path before trying to expand
#   Windows shortcuts.
# 2007-03-19
# o Now asterisk tags are handles dynamically, and not by the constructor.
# o Added getAsteriskTags().
# o Now the constructor expects a ChipEffectSetTuple instead of a list of
#   ChipEffectSet:s.
# o Updated code to internally make use of the ChipEffectSetTuple class.
# 2007-03-15
# o Updated the GladModel to only work with chromosome indices (integers).
# o Now the GladModel infers the set of possible chromosomes from the
#   GenomeInformation file.
# 2007-03-12
# o BUG FIX: getFullNames() of GladModel would give 'Error in getName(ceList[[
#   1]]) : no applicable method for "getName"' if there was not hybridization
#   for the *first* chip type in a set of multiple chip types.
# o BUG FIX: fit() of GladModel would give 'Error in FUN(X[[2]], ...) : no
#   applicable method for "getTags"' if there were not data files for all
#   chip types.
# 2007-02-20
# o Added getFullNames(), which for each tuple (across chip types) returns the
#   sample name of the tuple, together with all *common* tags across all
#   chip types.  Tags existing in only some of the chip types are ignored.
# 2007-02-16
# o Now the default version of the human genome is 'hg17' and not 'hg18'.
#   The reason for this is that the dChip annotation files are 'hg17'. We
#   still have to figure out a way to do version control for this in the
#   package.  Maybe it won't be a problem as soon as we start using the
#   annotation packages of Bioconductor.  On the to do list...
# o Added arguments 'organism' and 'db' to getRegions().
# 2007-02-15
# o Now getChipTypes() sorts the chip types in lexicographic order before
#   merging.  This guarantees the same result regardsless of order of the
#   input list.
# o Added getReportPath().
# o Path is now back to <rootPath>/<data set>,<tags>/<chipType>/.
# o Reports are written to reports/<data set>/<tags>/<chipType>/glad/.
# 2007-02-06
# o Updated the path to <rootPath>/<dataSetName>/<tags>/<chipType>/<set>/.
# 2007-01-25
# o Added so that plot() generates fixed sized horizontal margins (50px),
#   so that we can infer the relative genomic location from the horisontal
#   pixel location.
# 2007-01-21
# o Added a better error message when plot() fails to locate hgChromosomes.txt.
# 2007-01-17
# o Now argument 'arrays' can be either a vector of indices or array names,
#   or NULL.
# o Added indexOfArrays().
# 2007-01-16
# o Now NULL values for arguments 'arrays' and 'chromosomes' of fit() defaults
#   to all arrays and all chromosomes, respectively.
# o BUG FIX: writeRegions() would give an error if no regions was found.
# 2007-01-15
# o Now fit(..., force=TRUE) also calls getReferenceFiles(..., force=force).
# o Added some more Rdoc comments.
# 2007-01-10
# o Now plot() of GladModel is search for 'hgChromosomes.txt' in both
#   annotations/ and the package installation directory.
# 2007-01-07
# o Renamed MultiGladModel to GladModel fully replacing the older class.
# 2006-12-20
# o Now the class accepts any ChipEffectSet, not only CnChipEffectSet objects.
#   CnChipEffectSet objects are still validated specially, if used.
# 2006-12-17
# o BUG FIX: The new fitDone() in plot() choked on chr 23 (should be 'X').
# 2006-12-15
# o This class should be considered temporary, because we might design a
#   ChipEffectSet class that can contain multiple chip types, but treated as
#   if it contained one chip type, so it can be passed to the current 
#   GladModel class.  However, such a class design will require multiple 
#   inheritance etc, which will take time to develope.
# o Created from GladModel.R with history as below:
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
