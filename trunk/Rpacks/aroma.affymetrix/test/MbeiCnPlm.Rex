source("init.R")

verbose <- Arguments$getVerbose(TRUE);

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

units <- findUnitsTodo(model);
print(summary(units));

stop()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit the model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
units <- fit(model, units=56, force=TRUE, verbose=TRUE);
cat("Fitted ", length(units), " units.\n");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Calculate the average copy-number file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cesAvg <- getAverageFile(ces, verbose=TRUE);

ce <- getFile(ces, 1);
smoothScatterMvsA(ce, cesAvg, xlim=c(8,16));
abline(h=log(1:6/2, base=2), lty=3);
