library(oligo);
library(aroma.affymetrix);
source("init.R");
source("ui.R");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
choices <- c("MBEI,+-", "MBEI,A+B", "RMA,+-", "RMA,A+B");
ans <- menu(choices=choices, title="Select probe-level model (PLM)")
model <- choices[ans];
modelTags <- unlist(strsplit(model, split=","));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- c("normQuantile");
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
ds <- extract(dataSets[[1]], 1:7);
verbose && enter(verbose, "oligo");
verbose && print(verbose, ds);

verbose && enter(verbose, "Reading CEL files");
snpFSet <- extractSnpFeatureSet(ds, verbose=verbose);

verbose && cat(verbose, "Read data:");
verbose && print(verbose, snpFSet);
verbose && exit(verbose);

verbose && enter(verbose, "Preprocessing and summarizing array data");
snpQSet <- snprma(snpFSet);
verbose && cat(verbose, "Read data:");
verbose && print(verbose, snpQSet);
verbose && exit(verbose);

verbose && enter(verbose, "Fitting CRLMM");
fit <- crlmm(snpQSet);
verbose && exit(verbose);

verbose && exit(verbose);


