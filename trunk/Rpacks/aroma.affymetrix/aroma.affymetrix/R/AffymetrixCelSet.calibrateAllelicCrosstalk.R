###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod calibrateAllelicCrosstalk
#
# @title "Calibrates probepair signals for crosstalk between allele A and allele B"
#
# \description{
#  @get "title".  
# }
#
# @synopsis
#
# \arguments{
#   \item{outPath}{The path where to save the calibrated data files.}
#   \item{...}{Additional arguments passed to 
#     @seemethod "calibrateAllelicCrosstalk".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the calibrated @see "AffymetrixCelSet" object.
# }
#
# \details{
#  The calibration is done for each array separately.
# }
#
# \section{Benchmarking}{
#   On an IBM Thinkpad A31 with 1.8GHz and 1GB RAM:
#   \itemize{
#    \item{Mapping50K\_Xba240}{58960 SNPs: ~55-60s/array.}
#   }
# }
#
# @author
#
# \seealso{
#   @see "sfit::cfit".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("calibrateAllelicCrosstalk", "AffymetrixCelSet", function(this, path=NULL, name="calibAllelicCT", ..., verbose=FALSE) {
  require("sfit") || throw("Package 'sfit' not found.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  if (is.null(path)) {
    # Path structure: /calibAllelicCT/<data set name>/chip_data/<chip type>/
    path <- file.path(name, getName(this), "chip_data", getChipType(cdf));
  } 
  if (!is.null(path)) {
    # Verify this path (and create if missing)
    path <- Arguments$getWritablePath(path);
  }

  if (identical(getPath(this), path)) {
    throw("Cannot calibrate data file. Argument 'path' refers to the same path as the path of the data file to be calibrated: ", path);
  }
  mkdirs(path);


  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the cell indices for each possible allele basepair.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying cell indices for each possible allele basepair");
  setsOfProbes <- getAlleleProbePairs(cdf, verbose=verbose);
  gc <- gc();
  verbose && print(verbose, gc);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibrate each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calibrating ", nbrOfArrays(this), " arrays");
  verbose && enter(verbose, "Path: ", path);
  dataFiles <- list();
  nbrOfArrays <- nbrOfArrays(this);
  for (kk in seq_len(nbrOfArrays)) {
    df <- getFile(this, kk);
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                              kk, getName(df), nbrOfArrays));
    df <- calibrateAllelicCrosstalk(df, ..., path=path, 
                           setsOfProbes=setsOfProbes, verbose=less(verbose));
    dataFiles[[kk]] <- df;

    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  res <- newInstance(this, dataFiles);
  setCdf(res, cdf);

  res;
}, private=TRUE)


############################################################################
# HISTORY:
# 2006-10-06
# o make sure cdf association is inherited
# 2006-09-15
# o Adapted to the new aroma.affymetrix class structure.
# 2006-07-21
# o Added calibrateAllelicCrosstalk().
############################################################################
