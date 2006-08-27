source("init.R")

verbose <- Arguments$getVerbose(TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
path <- "chip_data/Hind/";
path <- "chip_data2/Xba/";
ds <- AffymetrixCelSet$fromFiles(path);

cdf <- getCdf(ds);

# Specify the model we want to fit
model <- AffymetrixCnRmaModel(ds, name="modelTotalCN");
print(model);

paf <- getProbeAffinities(model);
ces <- getChipEffects(model);
print(paf);
print(ces);

units <- findUnitsTodo(paf);
print(summary(units));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit the model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
units <- fit(model, moreUnits=5, verbose=TRUE);
cat("Fitted ", length(units), " units.\n");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Calculate the average copy-number file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cesAvg <- getAverageFile(ces, verbose=TRUE);

ce <- as.list(ces)[[1]];
smoothScatterMvsA(ce, cesAvg, xlim=c(8,16))
abline(h=log(1:6/2, base=2), lty=c(3,2,rep(3,4)))
