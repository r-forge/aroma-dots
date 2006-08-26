source("init.R")

verbose <- Arguments$getVerbose(TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
ds <- AffymetrixCelSet$fromFiles("cel/Hind/");

# Specify the model we want to fit
model <- AffymetrixLiWongModel(ds);
print(model);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit the model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit a few units.  The method is clever such that it will not 
# estimate the same units if not needed, instead it will read the
# data from the files where parameter estimates are stored.
f <- fit(model, units=50000:50010);
units <- fit(model, units=50000:50010, force=TRUE);
cat("Fitted ", length(units), " units.\n");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Investigated the fitted units
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Work with probe-affinity estimates only.  The object returned
# inherits from the AffymetrixCelFile class, that is, you can plot
# probe affinities spatially just as if they were probe signals 
# and so on.
paf <- getProbeAffinities(model);
print(paf);

