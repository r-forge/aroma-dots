for (zzz in 0) {

# Setup verbose output
verbose <- Arguments$getVerbose(-2)
timestampOn(verbose)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define an example dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Find any dataset
path <- findCelSet(minCount=5)
if (is.null(path))
  break

ds <- AffymetrixCelSet$fromFiles(path)
print(ds)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calibration
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsC <- calibrateAllelicCrosstalk(ds, verbose=verbose)
print(dsC)

} # for (zzz in 0)
rm(zzz)
