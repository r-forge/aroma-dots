###########################################################################/**
# @RdocClass CbsModel
#
# @title "The CbsModel class"
#
# \description{
#  @classhierarchy
#
#  This class represents the Circular Binary Segmentation (CBS) model [1]. 
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to the constructor of 
#              @see "CopyNumberSegmentationModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
#  @see "CopyNumberSegmentationModel".
# }
#
# \references{
#  [1] Olshen, A. B., Venkatraman, E. S., Lucito, R., Wigler, M. 
#      \emph{Circular binary segmentation for the analysis of array-based 
#      DNA copy number data. Biostatistics 5: 557-572, 2004.}\cr
#  [2] Venkatraman, E. S. & Olshen, A. B. 
#      \emph{A faster circular binary segmentation algorithm for the 
#      analysis of array CGH data}. Bioinformatics, 2007.\cr 
# }
#*/###########################################################################
setConstructorS3("CbsModel", function(cesTuple=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(cesTuple)) {
    require("DNAcopy") || throw("Package not loaded: DNAcopy");
  }

  extend(CopyNumberSegmentationModel(cesTuple=cesTuple, ...), "CbsModel")
})


setMethodS3("getAsteriskTag", "CbsModel", function(this, ...) {
  "CBS";
}, protected=TRUE)



setMethodS3("plot", "CbsModel", function(x, ..., pixelsPerMb=3, zooms=2^(0:7), pixelsPerTick=2.5, height=400, xmargin=c(50,50), imageFormat="current", skip=TRUE, path=NULL, callList=NULL, verbose=FALSE) {
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get genome annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for user specific annotation file
  filename <- "hgChromosomes.txt";
  pathname <- file.path("annotations", filename);
  if (!isFile(pathname)) {
    # If not found, fall back to the one in the package.
    pathname <- system.file("annotations", filename, 
                                                 package="aroma.affymetrix");
    if (!isFile(pathname))
      throw("Failed to locate file: ", filename);
  }

  genome <- readTable(pathname, header=TRUE, 
                            colClasses=c(nbrOfBases="integer"), row.names=1);

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
    setHook("onFit.fitGlad.CbsModel", NULL, action="replace");
  })

  setHook("onFit.fitGlad.CbsModel", function(fit, fullname) {
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
      cat("ERROR caught in onFit.fitGlad.CbsModel():\n");
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




##############################################################################
# HISTORY:
# 2007-08-20
# o Created from GladModel.R.
##############################################################################
