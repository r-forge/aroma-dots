for (zzz in 0) {

# Setup verbose output
verbose <- Arguments$getVerbose(-2)
timestampOn(verbose)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Find an example SNP-chip dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Search for a chip_data/<SNP chip type>/ directory...
for (dir in c("Hind", "Xba", "Nsp", "Sty")) {
  path <- filePath("chip_data", dir, expandLinks="any")
  if (isDirectory(path)) {
    # Try to define a dataset
    tryCatch({
      ds <- AffymetrixCelSet$fromFiles(path)
      verbose && cat(verbose, "Found data: ", path)
      break
    }, error = function(ex) {})
  }
  path <- NULL
}
if (is.null(path))
  break


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define an example dataset using this path
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("ds")) {
  ds <- AffymetrixCelSet$fromFiles(path)
}
verbose && ruler(verbose)
verbose && print(verbose, ds)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create a set of various PLMs for this dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("models", mode="list")) {
  mergeStrands <- TRUE
  combineAlleles <- TRUE
  models <- list(
    rma = RmaCnPlm(ds, mergeStrands=mergeStrands, 
                                      combineAlleles=combineAlleles),
    mbei = MbeiCnPlm(ds, mergeStrands=mergeStrands, 
                                      combineAlleles=combineAlleles),
  )
}
verbose && ruler(verbose)
verbose && print(verbose, models)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# For each model, fit a few units
#
# Note, by fitting the same set of units across models, the internal
# caching mechanisms of aroma.affymetrix makes sure that the data is
# only read into memory once. See log for reading speed.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
units <- 55+1:100

for (model in models) {
  verbose && ruler(verbose)
  unitsFitted <- fit(model, units=units, force=TRUE, verbose=verbose)
}
