###########################################################################/**
# @RdocClass MatSmoothing
#
# @title "The MatSmoothing class"
#
# \description{
#  @classhierarchy
#
#  This class represents a function for smoothing data with a trimmed mean.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeLevelTransform".}
#   \item{design}{A design @matrix.}
#   \item{probeWindow}{Bandwidth to use.  Effectively the width is 
#      2*probeWindow since it looks probeWindow bases in either direction.}
#   \item{nProbes}{The minimum number of probes to calculate a MAT score for.}
#   \item{meanTrim}{The amount of trimming of the mean in [0,0.5].}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \author{Mark Robinson (mrobinson[at]wehi.edu.au).}
#*/###########################################################################
setConstructorS3("MatSmoothing", function(..., design=NULL, probeWindow=300, nProbes=10, meanTrim=0.1) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(design)) {
    # Argument 'design':
    if (!is.matrix(design)) {
      throw("Argument 'design' is not a matrix: ", class(design)[1]);
    }

    # Argument 'probeWindow':
    probeWindow <- Arguments$getNumeric(probeWindow, range=c(0,Inf));

    # Argument 'nProbes':
    nProbes <- Arguments$getInteger(nProbes, range=c(1,Inf));

    # Argument 'meanTrim':
    meanTrim <- Arguments$getNumeric(meanTrim, range=c(0,0.5));
  }


  this <- extend(ProbeLevelTransform(...), "MatSmoothing",
    .design = design,
    .probeWindow = probeWindow,
    .nProbes = nProbes,
    .meanTrim = meanTrim
  );


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Further validation of the design matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getInputDataSet(this);
  if (!is.null(ds)) {
    # Validate the dimension of the design matrix
    nbrOfFiles <- nbrOfFiles(ds);
    design <- this$.design;
    dim <- dim(design);
    if (dim[1] != nbrOfFiles) {
      throw("The number of rows in the 'design' matrix, does not match the number of arrays in the input data set: ", dim[1], " != ", nbrOfFiles);
    }

#    if (dim[1] != dim[2]) {
#      throw("The 'design' matrix must be a square matrix: ", dim[1], "x", dim[2]);
#    }

    # Validate the contents of the design matrix
    for (cc in seq(length=ncol(design))) {
      if (!any(design[,cc] != 0)) {
        throw("Column #", cc, " in argument 'design' is all zero.");
      }
    } # for (cc ...)

    # Validate the names attributes of the design matrix
    outputNames <- colnames(design);
    if (is.null(outputNames)) {
      throw("Matrix 'design' does not have column names.");
    }
    outputNames <- Arguments$getCharacters(outputNames);
    if (any(duplicated(outputNames))) {
      throw("Argument 'design' contains duplicated column names: ", 
             paste(outputNames[duplicated(outputNames)]), collapse=", ");
    }

    # Check if the column names translates to valid filenames
    fullnames <- lapply(outputNames, FUN=function(fullname) {
      Arguments$getFilename(fullname);
    });
  }

  this;
})


setMethodS3("getAromaCellPositionFile", "MatSmoothing", function(this, ..., force=FALSE) {
  acp <- this$.acp;

  if (force || is.null(acp)) {
    dataSet <- getInputDataSet(this);
    cdf <- getCdf(dataSet);
    chipType <- getChipType(cdf, fullname=FALSE);
    nbrOfCells <- nbrOfCells(cdf);
    acp <- AromaCellPositionFile$byChipType(chipType, nbrOfCells=nbrOfCells, ...);
    this$.acp <- acp;
  }

  acp;
}, protected=TRUE)


setMethodS3("getParameters", "MatSmoothing", function(this, ...) {
  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, ...);

  # Get parameters of this class
  params2 <- list(
    design = this$.design,
    probeWindow = this$.probeWindow,
    nProbes = this$.nProbes,
    meanTrim = this$.meanTrim
  );

  # Append the two sets
  params <- c(params, params2);

  params;
}, private=TRUE)



setMethodS3("getExpectedOutputFiles", "MatSmoothing", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Match to the column names of the design matrix");
  params <- getParameters(this);
  design <- params$design;

  verbose && cat(verbose, "Expected result names:");
  fullnames <- colnames(design);
  verbose && str(verbose, fullnames);

  # Sanity check (backup)
  stopifnot(!is.null(fullnames));

  # "Dummy" filenames
  filenames <- sprintf("%s.CEL", fullnames);

  verbose && exit(verbose);

  fullnames;
}, protected=TRUE)




###########################################################################/**
# @RdocMethod process
#
# @title "Processes the data set"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, data already processed is re-processed,
#       otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "MatSmoothing", function(this, ..., units=NULL, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # calculates null distribution from first taking non-overlapping probes,
  # then replicating the negative half of the distribution.
  calcNullDist <- function(ch, ps, x) {
    MIN <- -999999
    #inds <- 
    y <- rep(MIN,length(x))
    n <- length(ch)
    indices <- split(seq_len(n), ch)
    nChr <- length(indices)
    count <- 0
    for (ii in seq_len(nChr)) {
      ind <- indices[[ii]]
      nInd <- length(ind)
      pos <- ind[1]
      for (jj in seq_len(nInd)) {
        pp <- ps[ ind[jj] ]
        if ( (pp-pos) > (probeWindow*2) ) {
          count <- count + 1
          y[count] <- x[ ind[jj] ]
          #inds[count] <- ind[jj]
          pos <- pp
        }  
      } # for (jj ...)
    } # for (ii ...)
    y <- y[y > MIN]
    md <- median(y)
    y <- y[y <= md]
    #list( m=md, sd = sd( c(y,-y+2*md) ), inds=inds[inds > MIN])
    list( m=md, sd = sd( c(y,-y+2*md) ) )
  } # callNullDist()
    
  # ------------------------------------------------------
  # tolerance function to be used below to get indices
  # ------------------------------------------------------
  tolFun <- function(u, tol) {
    which(abs(u) <= tol);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "MAT smoothing according to design matrix");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force && isDone(this)) {
    verbose && cat(verbose, "Already smoothed");
    verbose && exit(verbose);
    outputDataSet <- getOutputDataSet(this);
    return(invisible(outputDataSet));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get input data set
  ds <- getInputDataSet(this);

  # Get (and create) the output path
  outputPath <- getPath(this);
  
  # Get cdf
  cdf <- getCdf(ds);
  
  # Get which units to fit
  if (is.null(units)) {
    units <- seq_len(nbrOfUnits(cdf));
  }

  verbose && enter(verbose, "Locating probe position data");
  # Locate AromaCellPositionFile holding probe sequences
  acp <- getAromaCellPositionFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  # Get algorithm parameters
  params <- getParameters(this, verbose=less(verbose, 50));
  probeWindow <- params$probeWindow;
  design <- params$design;
  nProbes <- params$nProbes;
  meanTrim <- params$meanTrim;
  rm(params);
  
  # ------------------------------------------------------
  # internal function to do trimmed mean smoothing 
  # ------------------------------------------------------
  calcSmoothed <- function(posVector, dataMatrix, probeWindow, nProbes, meanTrim) {
    require(gsmoothr)
    nc <- ncol(dataMatrix)
    posM <- matrix(rep(posVector,nc), nc=nc)
    o <- order(posM)  # calculate ordering
	
    subsetInd <- seq(1,length(o), by=nc)
	smoothedScore <- tmeanC(posM[o], dataMatrix[o], probeWindow=probeWindow, nProbes=nProbes, trim=meanTrim)
	
	return(smoothedScore[subsetInd])
  }

  
  # ------------------------------------------------------
  # get CEL indices for units 
  # ------------------------------------------------------
  cdfIndices <- getCellIndices(cdf, stratifyBy="pm", units=units);

  # sanity check to validate that each group for tiling array has this 1 element
  nbrGroupsPerUnit <- nbrOfGroupsPerUnit(cdf)
  stopifnot(all(nbrGroupsPerUnit==1))
  
  cdfIndices <- lapply(cdfIndices, FUN=function(unit) unit$groups[[1]]$indices);
  unitNames <- names(cdfIndices);
  names(cdfIndices) <- NULL;

  nRows <- base::sapply(cdfIndices, FUN=length);
  allInds <- unlist(cdfIndices, use.names=FALSE);

  nbrOfUnits <- length(units);

  fullnamesOut <- colnames(design);
  # Sanity check (backup)
  stopifnot(!is.null(fullnamesOut));
  verbose && cat(verbose, "Result/output names:");
  verbose && str(verbose, fullnamesOut);

  
  # --------------------------------------------------------
  # loop through the number of columns in design matrix
  # calculate smoothed score for each combination of samples
  # --------------------------------------------------------
  for (ii in seq_len(ncol(design))) {
    fullname <- fullnamesOut[ii];
    verbose && enter(verbose, sprintf("Result file #%d ('%s') of %d", 
                                                  ii, fullname, ncol(design)));

    filename <- sprintf("%s.CEL", fullname);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

    # Already done?
    if (!force && isFile(pathname)) {
      verbose && cat(verbose, "Already processed.");
      verbose && exit(verbose);
      next;
    }

    # Check already here if there is already a tempory output file.  This
    # may be a left over from a interrupted previous run, or the fact that
    # the array is processed by another session elsewhere. /HB
    pathnameT <- sprintf("%s.tmp", pathname);
    pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=TRUE);

    matScoreNeg <- matScorePos <- nbrRows <- outputList <- vector("list", nbrOfUnits);
    names(outputList) <- unitNames;
  
    sampsKeep <- which(design[,ii] != 0);

    # Sanity check (already done in the setup checks, but a 2nd backup)
    stopifnot(length(sampsKeep) > 0);

    posOrNeg <- design[sampsKeep,ii];
      
    verbose && enter(verbose, "Reading probe data for the samples needed");
    dsII <- extract(ds, sampsKeep);
    dataList <- readUnits(dsII, units=units, stratifyBy="pm", verbose=verbose);
    rm(dsII);
    verbose && exit(verbose);
    
    dataList <- base::lapply(dataList, FUN=function(u) {
      matrix(log2(u[[1]]$intensities), ncol=length(sampsKeep))
    });
    
    verbose && enter(verbose, "Computing trimmed means for all units");
    # ------------------------------------------------------
    # loop through each unit
    # ------------------------------------------------------
    for (jj in seq_len( nbrOfUnits )) {
    
      # allocate space first time through
      if ( is.null( outputList[[jj]] ) ) {
        zeroes <- rep(0, times=nRows[jj]);
        matScorePos[[jj]] <- zeroes;
        matScoreNeg[[jj]] <- zeroes; 
        nbrRows[[jj]] <- zeroes;
        outputList[[jj]] <- zeroes;
      }
      if (jj %% 1000 == 0) {
        verbose && cat(verbose, sprintf("Completed %d/%d units ...", jj, nbrOfUnits));
      }
      
      indices <- cdfIndices[[jj]]
      
      # loop through all columns in matrix
      sampPos <- which(posOrNeg > 0)
      sampNeg <- which(posOrNeg < 0)
      nPos <- length(sampPos)
      nNeg <- length(sampNeg)
	  
      # extract probe positions
      pos <- acp[indices,2,drop=TRUE]
	  
      # if samples to smooth, do smoothing
      if( nPos )
	    matScorePos[[jj]] <- calcSmoothed(pos, dataList[[jj]][,sampPos,drop=FALSE], probeWindow=probeWindow, nProbes=nProbes, meanTrim=meanTrim)
      if( nNeg )
        matScoreNeg[[jj]] <- calcSmoothed(pos, dataList[[jj]][,sampNeg,drop=FALSE], probeWindow=probeWindow, nProbes=nProbes, meanTrim=meanTrim)
	}

    # Memory cleanup
    rm(dataList);

    verbose && exit(verbose);

    
    # calculate null distributions for scale factors
    if (length(sampNeg) > 0) {
      verbose && enter(verbose, "Gathering common info for calculating null distributions");
      chr <- acp[allInds,1, drop=TRUE];
      pos <- acp[allInds,2, drop=TRUE];
      verbose && exit(verbose);
      
      verbose && enter(verbose, "Calculating null distribution for controls");
      nullX <- unlist(matScoreNeg, use.names=TRUE);
      nullDistNeg <- calcNullDist(chr, pos, nullX);
      verbose && exit(verbose);
      
      verbose && enter(verbose, "Calculating null distribution for treatments");
      nullX <- unlist(matScorePos, use.names=TRUE);
      nullDistPos <- calcNullDist(chr, pos, nullX);
      verbose && exit(verbose);
      
      scaleFactor <- nullDistPos$sd / nullDistNeg$sd;
      
      # Memory cleanup
      rm(chr, nullX, nullDistNeg, nullDistPos);
      gc <- gc();
    } else {
      scaleFactor <- 1;
    } # if (length(sampNeg) > 0)

    
    for (jj in seq_len( nbrOfUnits )) {
      outputList[[jj]] <- matScorePos[[jj]] - scaleFactor*matScoreNeg[[jj]];
    }

    # Memory cleanup
    rm(matScoreNeg, matScorePos);
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Storing results
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Always validate the output file just before writing
    pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=TRUE);

    # Create CEL file to store results, if missing
    verbose && enter(verbose, "Creating CEL file for results, if missing");
    df <- getFile(ds, 1);
    createFrom(df, filename=pathnameT, path=NULL, verbose=less(verbose));
    verbose && exit(verbose);
    
    matScores <- unlist(outputList, use.names=FALSE);
    matScores <- 2^matScores;
    updateCel(pathnameT, indices=allInds, intensities=matScores, verbose=TRUE);

    # ...rename temporary file
    file.rename(pathnameT, pathname);
    if (!isFile(pathname) || isFile(pathnameT)) {
      throw("Failed to rename temporary file: ", pathnameT, " -> ", pathname);
    }
    rm(filename, pathname, pathnameT);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Next column in design matrix
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Memory cleanup
    rm(matScores, outputList, nbrRows);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # for (ii ...)

  # Clean up
  rm(cdfIndices, allInds);
  
  outputDataSet <- getOutputDataSet(this, force=TRUE);

  verbose && exit(verbose);
  
  invisible(outputDataSet);
})


############################################################################
# HISTORY:
# 2009-05-27 [MR]
# o 
# o Added sanity check for CDF file that all units have exactly 1 group
# o Added 'stratifyBy="pm"' to the call to readUnits() and getCellIndices()
#   since the smoothing is only to be done on PM probes
# 2009-05-25 [HB]
# o Added getExpectedOutputFiles() special to MatSmoothing. Removed 
#   isDone().
# 2009-05-23 [HB]
# o Added getOutputFiles() to MatSmoothing to work with new AromaTransform.
# o Now process() of MatSmoothing skips already process output files.
# o Updated process() of MatSmoothing to first write to a temporary file
#   which is then renamed.  This lower the risk for corrupt output files
#   due to processing interrupts.
# 2009-03-20 [MR]
# o Corrected the checking of isDone() for MatSmoothing objects.
# o If MatSmoothing has already been run, it returns the 
#   AffymetrixCelSet object.
# 2009-01-13 [HB]
# o MEMORY CLEANUP: Cleaning out more "done" variables and earlier.
# o Code cleanup.
# 2008-03-21
# o Created from BackgroundCorrection.R.
############################################################################
