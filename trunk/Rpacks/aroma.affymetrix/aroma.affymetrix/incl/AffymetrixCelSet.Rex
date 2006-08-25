
for (zzz in 0) {

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Look for the 100K Hind data directory
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- filePath("cel/Xba/", expandLinks="any")
if (!isDirectory(path))
  break


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define a dataset object for all files therein
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ds <- AffymetrixDataSet$fromFiles(path=path)

# Keep at most sex arrays for this example
ds <- extract(ds, 1:min(6,nbrOfArrays(ds)))
print(ds)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsQN <- normalizeQuantile(ds, subsetToAvg=0.50,
                                  typesToUpdate="pm", verbose=TRUE)
print(dsQN)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot densities
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
layout(matrix(1:2, ncol=1))
xlim <- c(6,16)
ylim <- c(0,0.40)
plotDensity(ds, types="pm", xlim=xlim, ylim=ylim, main="Raw data")
plotDensity(dsQN, types="pm", xlim=xlim, ylim=ylim, main="Normalized data")

# Clean up
layout(1)
rm(path, xlim, ylim)
# rm(ds, dsQN)
} # for (zzz in 0)
