source("init.R")

verbose <- Arguments$getVerbose(TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
path <- "chip_data2/Xba/";
ds <- AffymetrixCelSet$fromFiles(path);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create a set of models to work with
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("models", mode="list")) {
  models <- list(
    mbei = MbeiPlm(ds),
    rma = RmaPlm(ds),
    affine = AffinePlm(ds)
  )
}

mm <- 1;
model <- models[[mm]];
paf <- getProbeAffinities(model);
ces <- getChipEffects(model);
print(ces);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit the model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit a few units.  The method is clever such that it will not 
# estimate the same units if not needed, instead it will read the
# data from the files where parameter estimates are stored.
units <- fit(model, units=50000:50010, verbose=TRUE);
cat("Fitted ", length(units), " units.\n");

unit <- 50000;
print(paf[[unit]]);
print(ces[[unit]]);
