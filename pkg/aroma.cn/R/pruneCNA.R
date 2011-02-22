setMethodS3("pruneCNA", "PairedPSCBS", function(fit, ..., maxGeneration=Inf, onAtomicIsland=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'maxGeneration':
  maxGeneration <- Arguments$getDouble(maxGeneration, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
 

  fitT <- fit;

  fitList <- list();
  gg <- 1;
  while (gg <= maxGeneration) {
    verbose && enter(verbose, sprintf("Generation #%d", gg));

    fitList[[gg]] <- fitT;
  
    nbrOfSegments <- nbrOfSegments(fitT);
    maxH <- nbrOfSegments-2L;

    # Done
    if (maxH < 0) {
      break;
    }

    hasChanged <- FALSE;
    for (hh in 0:maxH) {
#    for (hh in 0) {
      verbose && enter(verbose, sprintf("Block size H=%d of %d", hh, maxH));
    
      res <- findAtomicAberrations(fitT, H=hh, ..., verbose=verbose);
    
      # (i) Atomic islands?
      if (hh == 0) {
        atomicIslands <- res$atomicRegions;
      } else {
        atomicIslands <- res$atomicIslands;
      }
      if (length(atomicIslands) > 0) {
        verbose && printf(verbose, "Atomic islands found (H=%d):\n", hh);
        verbose && print(verbose, atomicIslands);
    
        # Drop atomic islands and merge flanking segments
        dropList <- list();
        atomicIslands <- sort(atomicIslands, decreasing=TRUE);
        for (kk in seq(along=atomicIslands)) {
          atomicIsland <- atomicIslands[kk];
          verbose && enter(verbose, sprintf("Atomic island %d (%d) of %d", 
                           kk, atomicIsland, length(atomicIslands)));
          verbose && cat(verbose, "Number of segments before: ", nbrOfSegments(fitT));
          fitTT <- fitT;


          fitTT <- dropByRegions(fitTT, regions=atomicIsland, H=hh);

          fitDrop <- fitTT$dropped;
          fitTT$dropped <- NULL;
          dropList[[kk]] <- fitDrop;
    
          # Merge
          fitTT <- mergeTwoSegments(fitTT, left=atomicIsland-1L);

          if (is.function(onAtomicIsland)) {
            onAtomicIsland(fit0=fitT, fit1=fitTT, fitD=fitDrop,
                           atomicIsland=atomicIsland, H=hh);
          }

          fitT <- fitTT;

          verbose && cat(verbose, "Number of segments after: ", nbrOfSegments(fitT));
          verbose && exit(verbose);
        } # for (kk ...)
  
        fitT$dropped <- dropList;
        fitT$atomicIslands <- rev(atomicIslands);
        fitT$H <- hh;

        hasChanged <- TRUE;

        verbose && exit(verbose);
    
        # Go to next generation
        break;
      }
    
      verbose && exit(verbose);
    } # for (hh ...)
  
    if (!hasChanged) {
      break;
    }
  
    # Next generation
    gg <- gg + 1L;

    verbose && exit(verbose);
  } # while(gg <= maxGeneration)

  class(fitList) <- c("PruneCNA", class(fitList));

  fitList;
})



############################################################################
# HISTORY:
# 2011-01-18
# o Added class 'PruneCNA' to the return object of pruneCNA().
# o Now pruneCNA() returns pruned objects with the dropped segments 
#   included in a separate list.
# o Added argument 'maxGeneration' to pruneCNA().
# 2010-09-08
# o Updated to use findAtomicAberrations().
# 2010-07-24
# o CLEAN UP: Now the notation of the code better reflect the algorithm.
# o Now findAtomicRegions() returns ambigous atomic regions too.
# o Added argument 'ylim'.
# 2010-07-20
# o Added argument 'debugPlot'.
# 2010-07-19
# o Added trial version of segmentByPruneCBS().
# o TO DO: Down-weight loci that were close to earlier 
#   change points in the succeeding segmentations.
# o Added prototype version of findAtomicRegions().
# o Added prototype version of callByPruning().
# o Created.
############################################################################