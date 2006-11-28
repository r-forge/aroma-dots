library(aroma.affymetrix);
source("init.R");
source("ui.R");

pngDev <- System$findGraphicsDevice();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select what to do
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
choices <- c("Get GLAD model", "Fit GLAD", "Plot GLAD estimates", "Write GLAD table", "Write WIG tracks (for the UCSC Genome Browser)");
methods <- selectMenu(choices=choices, title="Select what you want to do (leave all blank to use all)");
if (length(methods) == 0) {
  methods <- seq(choices);
} else {
  methods <- match(methods, choices);
} 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- c("modelRmaCnPlm", "modelMbeiCnPlm", "normFragmentLength");
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
if (!exists("arrays"))
  arrays <- NULL;
arrays <- selectMenu(choices=choices, selected=arrays, title="Select the arrays you which to fit (leave all blank to use all)");
if (length(arrays) == 0) {
  arrays <- seq(choices);
} else {
  arrays <- match(arrays, choices);
} 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select chromosomes to fit
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
choices <- c(1:22,"X");
if (!exists("chromosomes"))
  chromosomes <- NULL;
chromosomes <- selectMenu(choices=choices, selected=chromosomes, title="Select the chromosomes you which to fit (leave all blank to use all)");
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


if (all(c(2,3) %in% methods) || all(c(2,4) %in% methods)) {
  methods <- setdiff(methods, 1);
}

res <- lapply(dataSets, FUN=function(ces) {
  verbose && enter(verbose, "Chip-effect set ", getName(ces));

  verbose && cat(verbose, "Chip-effect set:");
  verbose && print(verbose, ces);

  # Hard wired for now
  setMergeStrands(ces, TRUE);
  setCombineAlleles(ces, TRUE);

  # Define the GLAD model
  glad <- GladModel(ces);

  if (2 %in% methods) {
    verbose && enter(verbose, "Fitting GLAD");
    fit(glad, arrays=arrays, chromosomes=chromosomes, .retResults=FALSE, verbose=less(verbose));
    verbose && exit(verbose);
  }

  if (3 %in% methods) {
    verbose && enter(verbose, "Plotting GLAD (and fitting if not done)");
    plot(glad, arrays=arrays, chromosomes=chromosomes, verbose=less(verbose));
    verbose && exit(verbose);
  }

  if (4 %in% methods) {
    verbose && enter(verbose, "Writing GLAD regions to tabular file");
    # Filter out regions with only little change
    smoothing <- c(-Inf,-0.15,+0.15,+Inf);
    writeRegions(glad, arrays=arrays, chromosomes=chromosomes, smoothing=smoothing, oneFile=TRUE, skip=FALSE, verbose=less(verbose));
    verbose && exit(verbose);
  }

  if (5 %in% methods) {
    verbose && enter(verbose, "Writing WIG tracks");
    writeWig(glad, arrays=arrays, chromosomes=chromosomes, verbose=less(verbose));
    verbose && exit(verbose);
  }

  glad;
})
verbose && exit(verbose);

############################################################################
# HISTORY:
# 2006-11-xx
# o Recreated.
############################################################################
