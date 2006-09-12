source("init.R")

verbose <- Arguments$getVerbose(TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
path <- "chip_data3/Xba/";
ds <- AffymetrixCelSet$fromFiles(path);

# Specify the model we want to fit
model <- RmaSnpPlm(ds);
print(model);

ces <- getChipEffects(model);
paf <- getProbeAffinities(model);
ce <- getFile(ces, 1);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit the model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit a few units.  The method is clever such that it will not 
# estimate the same units if not needed, instead it will read the
# data from the files where parameter estimates are stored.
units <- fit(model, units=55+1:10, verbose=TRUE);
cat("Fitted ", length(units), " units.\n");

