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
#      2*probeWindow since it looks probeWindow bases in either direction}
#   \item{nProbes}{The minimum number of probes to calculate a MAT score for.}
#   \item{meanTrim}{The amount of trimming of the mean.}
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
  # Argument 'design':
  if (!is.null(design)) {
    if (!is.matrix(design)) {
      throw("Argument 'design' is not a matrix: ", class(design)[1]);
    }

    for (cc in seq(length=ncol(design))) {
      if (!any(design[,cc] != 0)) {
        throw("Column #", cc, " in argument 'design' is all zero.");
      }
    }
  }


  extend(ProbeLevelTransform(...), "MatSmoothing",
    .design = design,
    .probeWindow = probeWindow,
    .nProbes = nProbes,
    .meanTrim = meanTrim
  )
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


setMethodS3("isDone", "MatSmoothing", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Checking if data set is \"done\"");

  pathnames <- getOutputFiles(this);
  if (length(pathnames) == 0) {
    verbose && cat(verbose, "NOT done. No output files found.");
    return(FALSE);
  }
  
  tags <- paste(getTags(this),collapse=",")
  inputDs <- getInputDataSet(this);
  cdf <- getCdf(inputDs)
  ds <- AffymetrixCelSet$fromName(getName(inputDs),cdf=cdf,verbose=verbose,tags=tags)
  
  params <- getParameters(this);
  design <- params$design;

  if (length(pathnames) < ncol(design) ) {
    verbose && cat(verbose, "NOT done. Too few output files: ", 
                                   ncol(design), " < ", nbrOfFiles(ds));
    return(FALSE);
  }

  if (length(pathnames) > ncol(design) ) {
    throw("Too many output files found: ", 
                                  ncol(design), " > ", nbrOfFiles(ds));
  }

  if ( length( intersect(colnames(design), getNames(ds)) ) != nbrOfFiles(ds) ) {
    warning("Column names of the design matrix do not match the names of the dataset")
  }

  verbose && cat(verbose, "Done. All output files are there: ", ncol(design));
    
  verbose && exit(verbose);

  return(TRUE);
})



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
  tolFun <- function(u,tol) {
    which( abs(u) <= tol )
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, range=c(1,nbrOfUnits(cdf)));
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
  cdf <- getCdf(ds)
  
  # Get which units to fit
  if ( is.null(units) ) {
    units <- seq_len(nbrOfUnits(cdf))
  }

  verbose && enter(verbose, "Locating probe position data");
  # Locate AromaCellPositionFile holding probe sequences
  acp <- getAromaCellPositionFile(this, verbose=less(verbose, 5));
  verbose && exit(verbose);

  # Get algorithm parameters
  params <- getParameters(this, verbose=less(verbose, 50));
  probeWindow <- params$probeWindow
  design <- params$design
  nProbes <- params$nProbes
  meanTrim <- params$meanTrim
  rm(params);

  
  # ------------------------------------------------------
  # get CEL indices for units 
  # ------------------------------------------------------
  cdfIndices <- getCellIndices(cdf, units=units);

  # ASSUMPTION: the following is ok for tiling arrays since all have nGroups=1
  # MR /2008-2009
  cdfIndices <- lapply(cdfIndices, FUN=function(unit) unit$groups[[1]]$indices);
  unitNames <- names(cdfIndices);
  names(cdfIndices) <- NULL;

  nRows <- base::sapply(cdfIndices, FUN=length);
  allInds <- unlist(cdfIndices, use.names=FALSE);

  nUnits <- length( units )
  
  # ------------------------------------------------------
  # loop through the number of columns in design matrix
  # calculate smoothed score for each combination of samples
  # ------------------------------------------------------
  for (ii in seq_len(ncol(design))) {

    matScoreNeg <- matScorePos <- nbrRows <- outputList <- vector("list", nUnits)
    names(outputList) <- unitNames;
  
    sampsKeep <- which( design[,ii] != 0 )

    posOrNeg <- design[sampsKeep,ii]
      
    verbose && enter(verbose, "Reading probe data for the samples needed");
    dd <- readUnits(extract(ds,sampsKeep), units=units, verbose=verbose)
    verbose && exit(verbose)
    
    dd <- base::lapply(dd, FUN=function(u) {
      matrix(log2(u[[1]]$intensities), ncol=length(sampsKeep))
    })
    
    verbose && enter(verbose, "Computing trimmed means for all units");
    # ------------------------------------------------------
    # loop through each unit
    # ------------------------------------------------------
    for (jj in seq_len( nUnits )) {
    
      # allocate space first time through
      if ( is.null( outputList[[jj]] ) ) {
        matScorePos[[jj]] <- matScoreNeg[[jj]] <- nbrRows[[jj]] <- outputList[[jj]] <- rep(0, nRows[jj])
      }
      if (jj %% 1000==1)
        verbose && cat(verbose, sprintf("Completed %d/%d units ...", jj, nUnits));
      
      indices <- cdfIndices[[jj]]
      
      # loop through all columns in matrix
      sampPos <- which(posOrNeg > 0)
      sampNeg <- which(posOrNeg < 0)
      nPos <- length(sampPos)
      nNeg <- length(sampNeg)
      
      # calculate the pairwise distance matrix
      pos <- acp[indices,2,drop=TRUE]
      dist <- outer( pos, pos, FUN="-" )
      ppsRows <- base::apply(dist, MARGIN=1, FUN=tolFun, tol=probeWindow)
      
      # for every row in matrix
      ddJJ <- dd[[jj]];
      for (rw in seq_len( nRows[jj] )) {

        rows <- ppsRows[[rw]]
            
        nbrRows[[jj]][rw] <- length(rows)
                
        if ( length(rows) >= nProbes  ) {
        
          rootLength <- sqrt( round(length(rows)*(1-2*meanTrim)) )
          if ( nPos ) {
            vPos <- ddJJ[rows,sampPos]
            matScorePos[[jj]][rw] <- matScorePos[[jj]][rw] + mean(vPos, trim=meanTrim)*rootLength*sqrt(nPos)
          }
          if ( nNeg ) {
            vNeg <- ddJJ[rows,sampNeg]
            matScoreNeg[[jj]][rw] <- matScoreNeg[[jj]][rw] + mean(vNeg, trim=meanTrim)*rootLength*sqrt(nNeg)
          }
          
        }
      } # for (rw ...)
    } # for (jj ...)
    verbose && cat(verbose);

    # Memory cleanup
    rm(dd, ddJJ);

    verbose && exit(verbose);

    
    # calculate null distributions for scale factors
    if (length(sampNeg) > 0) {
      verbose && enter(verbose, "Gathering common info for calculating null distributions");
      chr <- acp[allInds,1,drop=TRUE]
      pos <- acp[allInds,2,drop=TRUE]
      verbose && exit(verbose)    
      
      verbose && enter(verbose, "Calculating null distribution for controls");
      nullX <- unlist(matScoreNeg, use.names=TRUE)
      nullDistNeg <- calcNullDist(chr, pos, nullX)
      verbose && exit(verbose)    
      
      verbose && enter(verbose, "Calculating null distribution for treatments");
      nullX <- unlist(matScorePos, use.names=TRUE)
      nullDistPos <- calcNullDist(chr, pos, nullX)
      verbose && exit(verbose)    
      
      scaleFactor <- nullDistPos$sd / nullDistNeg$sd
      
      # Memory cleanup
      rm(chr, nullX, nullDistNeg, nullDistPos);
      gc <- gc();
    } else {
      scaleFactor <- 1
    }
    
    for (jj in seq_len( nUnits )) {
      outputList[[jj]] <- matScorePos[[jj]] - scaleFactor*matScoreNeg[[jj]];
    }
    # Memory cleanup
    rm(matScoreNeg, matScorePos);
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Storing results
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    filename <- sprintf("%s.CEL", colnames(design)[ii]);
    pathname <- Arguments$getWritablePathname(filename, path=outputPath, ...);

    # Create CEL file to store results, if missing
    verbose && enter(verbose, "Creating CEL file for results, if missing");
    df <- getFile(ds,1)
    createFrom(df, filename=pathname, path=NULL, verbose=less(verbose));
    verbose && exit(verbose);
    
    matScores <- unlist(outputList, use.names=FALSE);
    matScores <- 2^matScores;
    updateCel(pathname, indices=allInds, intensities=matScores, verbose=TRUE);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Next column in design matrix
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Memory cleanup
    rm(matScores, outputList, nbrRows);

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);
  } # for (ii ...)

  # Clean up
  rm(cdfIndices, allInds);
  
  outputDataSet <- getOutputDataSet(this, force=TRUE);

  verbose && exit(verbose);
  
  invisible(outputDataSet);
})


############################################################################
# HISTORY:
# 2009-03-20 [MR]
# o Corrected the checking of isDone() for MatSmoothing objects
# o If MatSmoothing has already been run, it returns the AffymetrixCelSet object
# 2009-01-13 [HB]
# o MEMORY CLEANUP: Cleaning out more "done" variables and earlier.
# o Code cleanup.
# 2008-03-21
# o Created from BackgroundCorrection.R.
############################################################################
