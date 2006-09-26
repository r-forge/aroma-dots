# Usage:
#  R --no-save --args --chiptype="HuEx-1_0-st-v2" < rConvertCdf.R

source("init.R")
rConvertCdf <- function(chipType=NULL, destpath="cdf", overwrite=FALSE, verbose=TRUE, ...) {
  require(R.utils) || throw("Package R.utils not loaded.");
  require(aroma.affymetrix) || throw("Package aroma.affymetrix not loaded.");

  args <- commandArgs(asValues=TRUE, excludeReserved=TRUE);
  print(args);

  if (!is.null(args$chiptype))
    chiptype <- args$chiptype;
  if (!is.null(args$destpath))
    destpath <- args$destpath;
  if (!is.null(args$overwrite))
    overwrite <- args$overwrite;
  if (!is.null(args$verbose))
    verbose <- args$verbose;

  # Argument 'chiptype':
  if (is.null(chiptype))
    throw("Argument 'chiptype' must not be empty: NULL");

  # Argument 'destpath':
  destpath <- Arguments$getWritablePath(destpath);

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite);

  filename <- sprintf("%s.cdf", chiptype);
  pathname <- Arguments$getWritablePathname(filename, path=destpath, 
                                                 mustNotExists=!overwrite);
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  timestampOn(verbose);

  verbose && enter(verbose, "Trying to convert '", chiptype, "' CDF");

  # Define the source CDF
  verbose && enter(verbose, "Defining CDF");
  cdf <- AffymetrixCdfFile$fromChipType(chiptype, verbose=verbose);
  print(cdf);
  verbose && exit(verbose);

  # Convert 
  verbose && enter(verbose, "Converting CDF");
  cdf2 <- convert(cdf, verbose=verbose)
  print(cdf);
  verbose && exit(verbose);

  verbose && exit(verbose);
} # rConvertCdf()


# Call the main function
rConvertCdf();

###########################################################################
# BENCHMARKING:
#
# 20060925 19:30:59|Trying to convert 'HuEx-1_0-st-v2.text' CDF...
# 20060925 19:30:59| Defining CDF...
# AffymetrixCdfFile:
# Path: /accounts/gen/vis/hb/data/Affymetrix/cdf
# Filename: HuEx-1_0-st-v2.text.cdf
# Filesize: 933.84Mb
# Chip type: HuEx-1_0-st-v2
# Dimension: 2560x2560
# Number of cells: 6553600
# Number of units: 1432154
# Cells per unit: 4.58
# Number of QC units: 0
# RAM: 0.00Mb
# 20060925 19:30:59| Defining CDF...done
# 20060925 19:30:59| Converting CDF...
# Reading CDF header...
# Reading CDF header...done
# Reading CDF QC units...
# Reading CDF QC units...done
# Reading CDF units...
# Reading CDF units...done
# Writing CDF structure...
# 20060925 19:39:49| Time: Mon Sep 25 19:39:49 PDT 2006
# [top: 4.5g (29.8%) RAM, 99.9% CPU]

#
# HISTORY:
# 2006-09-24
# o Created.
###########################################################################
