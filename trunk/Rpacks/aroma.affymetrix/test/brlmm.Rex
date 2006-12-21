library(oligo);
library(aroma.affymetrix);
source("init.R");
source("ui.R");

moreUnits <- 30;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- "probeData";
pattern <- "^[a-zA-Z0-9_-]*$";
pattern <- NULL;
paths <- "raw";
pattern <- "^[^,]*$";
pattern <- "^[a-zA-Z0-9_-]*$";
if (exists("dataSets", mode="list")) {
  ds <- dataSets[[1]];
  name <- paste(c(getName(ds), getTags(ds)), collapse=",");
  ans <- menu(c("<new data sets>", name), title="Select data set(s) to process:");
  rm(ds);
  if (ans == 1)
    dataSets <- selectDataSets(paths, pattern=pattern);
} else {
  dataSets <- selectDataSets(paths, pattern=pattern);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
verbose && enter(verbose, "Processing ", length(dataSets), " data sets");

res <- list()
for (kk in seq(along=dataSets)) {
  ds <- dataSets[[kk]];
  verbose && enter(verbose, sprintf("Data set %d", kk));
  verbose && print(verbose, ds);

  apt <- AptProbesetGenotype(ds);
  verbose && print(verbose, apt);

  callSet <- process(apt, verbose=verbose);
  verbose && print(verbose, callSet);

  res[[kk]] <- callSet;
}

verbose && exit(verbose);

