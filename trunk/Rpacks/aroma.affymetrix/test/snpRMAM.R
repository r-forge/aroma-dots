source("init.R")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdfStrandJoiner <- function(units) {
  if (is.list(units[[1]])) {
    applyCdfGroups(units, cdfMergeStrands);
  } else {
    lapply(units, FUN=unique);
  }
}

verbose <- Arguments$getVerbose(TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
ds <- AffymetrixCelSet$fromFiles("chip_data/Hind/");

# Tweak the CDF such that forward and reverse strands are treated
# together.  Thus, each time the CDF is queried the below function
# restructure the CDF layout.  An alternative is to create a custom
# designed CDF file and use that instead.
cdf <- getCdf(ds);
setRestructor(cdf, cdfStrandJoiner);


# Specify the model we want to fit
model <- AffymetrixRmaModel(ds);
print(model);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit the model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Fit a few units.  The method is clever such that it will not 
# estimate the same units if not needed, instead it will read the
# data from the files where parameter estimates are stored.
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

