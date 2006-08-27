source("init.R")

verbose <- Arguments$getVerbose(TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
path <- "chip_data/Hind/";
path <- "chip_data1/Nsp/";
ds <- AffymetrixCelSet$fromFiles(path);

cdf <- getCdf(ds);

# Specify the model we want to fit
model <- AffymetrixCnRmaModel(ds, name="modelTotalCN");
print(model);

paf <- getProbeAffinities(model);
ces <- getChipEffects(model);
print(paf);
print(ces);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit the model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
units <- fit(model, verbose=TRUE);
cat("Fitted ", length(units), " units.\n");


