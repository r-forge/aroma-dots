library(aroma.affymetrix);
source("init.R");
source("ui.R");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
pngDev <- System$findGraphicsDevice();
ext <- "png";
if (ext == "png") {
  height <- 640;
  dev <- pngDev;
} else if (ext == "eps") {
  height <- 8 / 1.618;
  dev <- eps;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- c("modelRmaCnPlm", "modelMbeiCnPlm");
pattern <- ",A[+]B";
if (exists("dataSets", mode="list") && inherits(dataSets[[1]], "CnChipEffectSet")) {
  ds <- dataSets[[1]];
  name <- paste(getName(ds), paste(getTags(ds), collapse=","), sep=",");
  ans <- menu(c("<new data sets>", name), title="Select data set(s) to process:");
  if (ans == 1)
    dataSets <- selectDataSets(paths, pattern=pattern, class=CnChipEffectSet);
} else {
  dataSets <- selectDataSets(paths, pattern=pattern, class=CnChipEffectSet);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select arrays 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ces <- dataSets[[1]];
choices <- getNames(ces);
if (!exists("arrays"))
  arrays <- NULL;
arrays <- selectMenu(choices=choices, selected=arrays, title="Select the arrays to normalize (leave all blank to use all)");
if (length(arrays) == 0) {
  arrays <- seq(choices);
} else {
  arrays <- match(arrays, choices);
} 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
verbose && enter(verbose, "Processing ", length(dataSets), " data sets");
s <- paste("Arrays: ", seqToHumanReadable(arrays), sep="");
s <- GenericSummary(s);
print(s);

dataSets <- lapply(dataSets, function(ces) {
  extract(ces, arrays);
})

res <- lapply(dataSets, function(ces) {
  verbose && enter(verbose, "Chip-effect set ", getName(ces));

  verbose && cat(verbose, "Chip-effect set:");
  verbose && print(verbose, ces);

  # Hard wired for now
  setMergeStrands(ces, TRUE);
  setCombineAlleles(ces, TRUE);

  norm <- FragmentLengthNormalization(ces, subsetToFit=1/3);
  verbose && print(verbose, norm);

  cesFLN <- process(norm, verbose=verbose);

  verbose && exit(verbose);

  cesFLN;
})
verbose && exit(verbose);


############################################################################
# HISTORY:
# 2006-11-28
# o Recreated.
############################################################################
