for (zzz in 0) {

# Setup verbose output
verbose <- Arguments$getVerbose(-3)
timestampOn(verbose)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Find an example dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("ds")) {
  # Search for a chip_data/ directory...
  dataPath <- filePath("chip_data", expandLinks="any")
  dirs <- list.files(path=dataPath, full.names=TRUE)
  ds <- NULL
  for (path in dirs) {
    if (isDirectory(path)) {
      # Try to define a dataset
      tryCatch({
        ds <- AffymetrixCelSet$fromFiles(path)
        verbose && cat(verbose, "Found a dataset: ", path)
        break
      }, error = function(ex) {})
    }
  }
  if (is.null(ds))
    break
}

verbose && cat(verbose, "Working with dataset:")
print(ds)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsQN <- normalizeQuantiles(ds, subsetToAvg=0.50, typesToUpdate="pm", verbose=verbose)
print(dsQN)


} # for (zzz in 0)
rm(zzz)
