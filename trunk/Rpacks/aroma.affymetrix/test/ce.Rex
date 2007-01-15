library(aroma.affymetrix);
source("init.R");
source("ui.R");

pngDev <- System$findGraphicsDevice();
skip <- TRUE;
force <- FALSE;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- "plmData";
pattern <- ",A[+]B";
if (exists("dataSets", mode="list") && inherits(dataSets[[1]], "CnChipEffectSet")) {
  ds <- dataSets[[1]];
  name <- paste(getName(ds), paste(getTags(ds), collapse=","), sep=",");
  ans <- menu(c("<new data sets>", name), title="Select data set(s) to process:");
  rm(ds);
  if (ans == 1)
    dataSets <- selectDataSets(paths, pattern=pattern, class=CnChipEffectSet);
} else {
  dataSets <- selectDataSets(paths, pattern=pattern, class=CnChipEffectSet);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Define GLAD model
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
glad <- GladModel(dataSets);
print(glad);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select arrays to fit
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
choices <- getNames(glad);
if (!exists("arrays", mode="integer"))
  arrays <- NULL;
arrays <- selectMenu(choices=choices, selected=arrays, title="Select the arrays you which to fit (leave all blank to use all)");
if (length(arrays) == 0) {
  arrays <- seq(along=choices);
} else {
  arrays <- match(arrays, choices);
} 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select chromosomes to fit
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
choices <- getChromosomes(glad);
if (!exists("chromosomes"))
  chromosomes <- NULL;
chromosomes <- selectMenu(choices=choices, selected=chromosomes, title="Select the chromosomes you which to fit (leave all blank to use all)");
if (length(chromosomes) == 0) {
  chromosomes <- choices;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Process
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
verbose && enter(verbose, "Generate ChromosomeExplorer report");
ce <- ChromosomeExplorer(glad);
print(ce);

s <- paste("Arrays: ", seqToHumanReadable(arrays), sep="");
nas <- is.na(as.integer(chromosomes));
str <- paste(seqToHumanReadable(chromosomes[!nas]), 
             paste(chromosomes[nas], collapse=","), sep=",");
s <- c(s, paste("Chromosomes: ", str, sep=""));
s <- GenericSummary(s);
print(s);

process(ce, arrays=arrays, chromosomes=chromosomes, skip=skip, 
                                                     verbose=less(verbose));
verbose && exit(verbose);

############################################################################
# HISTORY:
# 2007-01-04
# o Created.
############################################################################
