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
#   \item{setsOfProbes}{A named @list where each item specifies the probe 
#     indices for each of the four groups; forward and reverse strands 
#     on allele A and allele B.}
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
#    \item{Mapping50K_Xba240}{58960 SNPs: ~55-60s/array.}
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
setMethodS3("calibrateAllelicCrosstalk", "AffymetrixCelSet", function(this, path=NULL, name="calibAllelicCrosstalk", ..., setsOfProbes=NULL, verbose=FALSE) {
  require(sfit) || throw("Package 'sfit' not found.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);

  # Argument 'path':
  if (is.null(path)) {
    # Path structure: <data-set name>/<normalization name>/<chip type>/
    # Compare with  : <data-set name>/chip_data/<chip type>/
    path <- filePath(getName(this), name, getChipType(cdf));
  } 
  if (!is.null(path)) {
    path <- Arguments$getWritablePath(path);
  }
  mkdirs(path);
  if (identical(getPath(this), path)) {
    throw("Cannot calibrate data file. Argument 'path' refers to the same path as the path of the data file to be calibrated: ", path);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the cell indices for each possible allele basepair.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(setsOfProbes)) {
    verbose && enter(verbose, "Identifying cell indices for each possible allele basepair");
    setsOfProbes <- getAlleleProbePairs(cdf, verbose=verbose);
    gc();
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calibrate each array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calibrating ", nbrOfArrays(this), " arrays");
  verbose && enter(verbose, "Path: ", path);
  dataFiles <- list();
  for (kk in seq(this)) {
    df <- getFile(this, kk);
    verbose && enter(verbose, "Array #", kk, " (", getName(df), ")");
    df <- calibrateAllelicCrosstalk(df, ..., path=path, 
                                 setsOfProbes=setsOfProbes, verbose=less(verbose));
    dataFiles[[kk]] <- df;

    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  newInstance(this, dataFiles);
}) # calibrateAllelicCrosstalk()


############################################################################
# HISTORY:
# 2006-09-15
# o Adapted to the new aroma.affymetrix class structure.
# 2006-07-21
# o Added calibrateAllelicCrosstalk().
############################################################################
