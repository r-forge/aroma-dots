library(aroma.affymetrix)
source("init.R");
source("ui.R");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- "raw";
pattern <- "_(199|200)[0-9]";
if (exists("dataSets", mode="list")) {
  ds <- dataSets[[1]];
  name <- paste(getName(ds), paste(getTags(ds), collapse=","), sep=",");
  ans <- menu(c("<new data sets>", name), title="Select data set(s) to process:");
  if (ans == 1)
    dataSets <- selectDataSets(paths, pattern=pattern);
} else {
  dataSets <- selectDataSets(paths, pattern=pattern);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
verbose && enter(verbose, "Processing ", length(dataSets), " data sets");
dataSets <- lapply(dataSets, FUN=function(ds) {
  verbose && enter(verbose, "Data set ", getName(ds));

  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  calib <- AllelicCrosstalkCalibrator(ds);
  dsOut <- process(calib, verbose=less(verbose));

  verbose && cat(verbose, "Output data set:");
  verbose && print(verbose, dsOut);

  dsOut;
})
verbose && exit(verbose);

