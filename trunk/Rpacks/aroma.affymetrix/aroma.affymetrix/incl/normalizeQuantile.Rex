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
# Normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsQN <- normalizeQuantiles(ds, subsetToAvg=0.50,
                               typesToUpdate="pm", verbose=verbose)
print(dsQN)

} # for (zzz in 0)
rm(zzz)
