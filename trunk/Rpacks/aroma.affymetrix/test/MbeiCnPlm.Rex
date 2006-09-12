source("init.R")

verbose <- Arguments$getVerbose(TRUE);
timestampOn(verbose);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
path <- "chip_data3/Xba/";
ds <- AffymetrixCelSet$fromFiles(path);

# Specify the model we want to fit
model <- MbeiCnPlm(ds, mergeStrands=TRUE, combineAlleles=TRUE);
print(model);

paf <- getProbeAffinities(model);
ces <- getChipEffects(model);
print(paf);
print(ces);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit the model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
units <- fit(model, moreUnits=5, verbose=verbose);
cat("Fitted ", length(units), " units.\n");


stop("Done!")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Calculate the average copy-number file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cesAvg <- getAverageFile(ces, verbose=verbose);

ce <- getFile(ces, 1);
smoothScatterMvsA(ce, cesAvg, xlim=c(8,16));
abline(h=log(1:6/2, base=2), lty=3);
