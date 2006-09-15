# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Define a dataset object based on all CEL files in a directory
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "MyDataSet/chip_data/Xba/"
ds <- AffymetrixCelSet$fromFiles(path)
print(ds)

# Keep at most sex arrays for this example
ds <- extract(ds, 1:min(6,nbrOfArrays(ds)))
print(ds)

