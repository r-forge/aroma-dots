###########################################################################/**
# @RdocClass AlleleSummation
#
# @title "The AlleleSummation class"
#
# \description{
#  @classhierarchy
#
#  This class takes allele-specific chip effect estimates of a 
#  SnpChipEffectSet and returns a CnChipEffectSet holding the summed
#  allele estimates.
# }
# 
# @synopsis
#
# \arguments{
#   \item{dataSet}{A @see "SnpChipEffectSet".}
#   \item{...}{Arguments passed to @see "UnitModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("AlleleSummation", function(dataSet=NULL, ignoreNAs=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "SnpChipEffectSet")) {
      throw("Argument 'dataSet' is not a SnpChipEffectSet object: ", class(dataSet)[1]);

    }
  }

  extend(UnitModel(dataSet=dataSet, ...), "AlleleSummation",
    ignoreNAs = ignoreNAs,
    "cached:.outputSet" = NULL
  )
})


setMethodS3("getAsteriskTags", "AlleleSummation", function(this, collapse=NULL, ...) {
  # Returns 'U' (but allow for future extensions)
  tags <- NextMethod("getAsteriskTags", this, collapse=NULL);
  tags[1] <- "SA";

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
})


setMethodS3("getRootPath", "AlleleSummation", function(this, ...) {
  "plmData";
}, private=TRUE)



###########################################################################/**
# @RdocMethod getChipEffectSet
#
# @title "Gets the set of chip effects for this model"
#
# \description{
#  @get "title".
#  There is one chip-effect file per array.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "ChipEffectSet" object.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChipEffectSet", "AlleleSummation", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # The output set
  outputSet <- this$.outputSet;
  if (!is.null(outputSet))
    return(outputSet);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create chip-effect files 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Let the parameter object know about the CDF structure, because we 
  # might use a modified version of the one in the CEL header.
  ds <- getDataSet(this);
  if (length(ds) == 0)
    throw("Cannot create chip-effect set. The input data set is empty.");
  
  verbose && enter(verbose, "Getting chip-effect set from data set");
  cdfM <- getCdf(ds);

  # Gets the ChipEffects Class object
  clazz <- getChipEffectSetClass(this);
  outputSet <- clazz$fromDataSet(dataSet=ds, path=getPath(this), cdf=cdfM,
                                                    verbose=less(verbose));
  setMergeStrands(outputSet, getMergeStrands(ds));
  setCombineAlleles(outputSet, TRUE);
  verbose && exit(verbose);

  # Store in cache
  this$.outputSet <- outputSet;

  outputSet;
})


setMethodS3("getChipEffectSetClass", "AlleleSummation", function(static, ...) {
  CnChipEffectSet;
}, static=TRUE, private=TRUE)


setMethodS3("process", "AlleleSummation", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Summing allele-specific estimates");
  inputSet <- getDataSet(this);
  verbose && print(verbose, inputSet);

  outputSet <- getChipEffectSet(this);

  cdf <- getCdf(this);
  snps <- indexOf(cdf, "SNP");
  otherUnits <- setdiff(seq(length=nbrOfUnits(cdf)), snps);

  snpUgcMap <- otherUgcMap <- NULL;

  ignoreNAs <- this$ignoreNAs;

  nbrOfArrays <- nbrOfArrays(inputSet);
  for (aa in seq(length=nbrOfArrays)) {
    verbose && enter(verbose, sprintf("Array #%d of %d", aa, nbrOfArrays));
    inputFile <- getFile(inputSet, aa);
    outputFile <- getFile(outputSet, aa);

    verbose && enter(verbose, "Copying signals for non-SNP units");
    if (is.null(otherUgcMap)) {
      # getUnitGroupCellMap <- getCellMap;  # To do!
      otherUgcMap <- getCellMap(inputFile, units=otherUnits);
    }
    cells <- otherUgcMap$cell;
    data <- readCel(getPathname(inputFile), indices=cells, 
                    readIntensities=TRUE, readStdvs=TRUE, readPixels=TRUE);
    data <- as.data.frame(data[c("intensities", "stdvs", "pixels")]);
    verbose && str(verbose, data);
    data <- cbind(cell=cells, data);
    rm(cells);
    updateDataFlat(outputFile, data=data);
    rm(data);
    verbose && exit(verbose);

    verbose && enter(verbose, "Summing allele signals for SNP units");
    if (is.null(snpUgcMap)) {
      snpUgcMap <- getCellMap(inputFile, units=snps);
    }
    cells <- snpUgcMap$cell;
    data <- readCel(getPathname(inputFile), indices=cells, 
                    readIntensities=TRUE, readStdvs=FALSE, readPixels=FALSE);
    yAB <- data[["intensities"]];
    verbose && cat(verbose, "(yA,yB) signals:");
    verbose && str(verbose, yAB);
    # Next we will assume that the data points are ordered as allele
    # (A,B,A,B,A,B,...)
    yAB <- matrix(yAB, nrow=2);

    # Sum the alleles
    y <- rep(NA, ncol(yAB));
    okAB <- !is.na(yAB);
    # 1) No missing data
    ok <- okAB[1,] & okAB[2,];
    y[ok] <- yAB[1,ok] + yAB[2,ok];
    if (ignoreNAs) {
      # 2a) Missing data in allele A
      ok <- !okAB[1,] & okAB[2,];
      y[ok] <- yAB[2,ok];
      # 2b) Missing data in allele B
      ok <- okAB[1,] & !okAB[2,];
      y[ok] <- yAB[1,ok];
    }
    verbose && cat(verbose, "y=yA+yB signals:");
    verbose && str(verbose, y);
    # Store signals in the cell for the A alleles:
    cells <- matrix(cells, nrow=2);
    cells <- cells[1,];

    data <- cbind(cell=cells, intensities=y);
    rm(cells, y);

    updateDataFlat(outputFile, data=data);
    rm(data);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (aa ...)

  verbose && exit(verbose);

  outputSet;
})






############################################################################
# HISTORY:
# 2008-02-20
# o Created.
############################################################################
