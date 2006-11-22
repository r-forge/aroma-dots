library(aroma.affymetrix);
source("init.R");
source("ui.R");

pngDev <- System$findGraphicsDevice();

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
# Select arrays to fit
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ces <- dataSets[[1]];
choices <- getNames(ces);
arrays <- selectMenu(choices=choices, title="Select the arrays you which to fit (leave all blank to use all)");
if (length(arrays) == 0) {
  arrays <- seq(ces);
} else {
  arrays <- match(arrays, choices);
} 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select chromosomes to fit
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
choices <- c(1:22,"X");
chromosomes <- selectMenu(choices=choices, title="Select the chromosomes you which to fit (leave all blank to use all)");
if (length(chromosomes) == 0) {
  chromosomes <- choices;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
verbose && enter(verbose, "Processing ", length(dataSets), " data sets");
s <- paste("Arrays: ", seqToHumanReadable(arrays), sep="");
nas <- is.na(as.integer(chromosomes));
str <- paste(seqToHumanReadable(chromosomes[!nas]), 
             paste(chromosomes[nas], collapse=","), sep=",");
s <- c(s, paste("Chromosomes: ", str, sep=""));
s <- GenericSummary(s);
print(s);

res <- lapply(dataSets, FUN=function(ces) {
  verbose && enter(verbose, "Chip-effect set ", getName(ces));

  verbose && cat(verbose, "Chip-effect set:");
  verbose && print(verbose, ces);

  # Hard wired for now
  setMergeStrands(ces, TRUE);
  setCombineAlleles(ces, TRUE);

  # Define the GLAD model
  glad <- GladModel(ces);

  verbose && enter(verbose, "Fitting & plotting GLAD");
  plot(glad, arrays=arrays, chromosomes=chromosomes, verbose=less(verbose));
  verbose && exit(verbose);

  verbose && enter(verbose, "Writing CN regions to file");
  # Filter out regions with only little change
  smoothing <- c(-Inf,-0.15,+0.15,+Inf);
  writeRegions(glad, arrays=arrays, chromosomes=chromosomes, smoothing=smoothing, verbose=less(verbose));
  verbose && exit(verbose);

  glad;
})
verbose && exit(verbose);
