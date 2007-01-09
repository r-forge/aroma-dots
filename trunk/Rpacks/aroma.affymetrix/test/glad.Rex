library(aroma.affymetrix);
source("init.R");
source("ui.R");

pngDev <- System$findGraphicsDevice();
skip <- TRUE;
imageFormat <- "screen";
imageFormat <- "png";

force <- FALSE;
force <- TRUE;
suppressWarnings({rm(res)});

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select what to do
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
choices <- c("Get GLAD model", "Fit GLAD", "List SD", "Plot GLAD estimates", "Write GLAD table", "Write WIG tracks (for the UCSC Genome Browser)");
if (!exists("methods", mode="integer"))
  methods <- NULL;
methods <- selectMenu(choices=choices, selected=methods, title="Select what you want to do (leave all blank to use all)");
if (length(methods) == 0) {
  methods <- seq(choices);
  names(methods) <- choices;
} else {
  methods <- match(methods, choices);
  names(methods) <- choices[methods];
} 

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
for (kk in seq(along=dataSets)) {
  ces <- dataSets[[kk]];

  # Hard wired for now
  setMergeStrands(ces, TRUE);
  setCombineAlleles(ces, TRUE);

  dataSets[[kk]] <- ces;
  rm(ces);
}

glad <- GladModel(dataSets);


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
verbose && enter(verbose, "Fit GLAD to chip effects");
print(glad);

s <- paste("Chip types: ", paste(getChipTypes(glad), collapse=", "), sep="");
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

if (2 %in% methods) {
  verbose && enter(verbose, "Fitting GLAD");
  fit(glad, arrays=arrays, chromosomes=chromosomes, force=force, .retResults=FALSE, verbose=less(verbose));
  verbose && exit(verbose);
}

if (3 %in% methods) {
  verbose && enter(verbose, "Display SD");
  sds <- list();
  for (chr in chromosomes) {
    fit <- fit(glad, arrays=arrays, chromosomes=chr, force=force, .retResults=FALSE);
    sdEst <- lapply(fit, function(cc) cc[[1]]$SigmaC$Value);
    sdEst <- unlist(sdEst);
    sds[[chr]] <- sdEst;
  }
  print(sds);
  verbose && exit(verbose);
}

if (4 %in% methods) {
  verbose && enter(verbose, "Plotting GLAD (and fitting if not done)");

  # Check for genotype calls
  callList <- list();
  cesList <- getListOfChipEffects(glad);
  for (kk in seq(along=cesList)) {
    ces <- cesList[[kk]];
    path <- getPath(ces);
    path <- gsub("plm", "crlmm", path, fixed=TRUE);
    path <- gsub(",QN,RMA,A+B,FLN", ",PPC,OQN,RMA,oligo,+-", path, fixed=TRUE);
    if (isDirectory(path)) {
      verbose && cat(verbose, "Detected genotype calls");
      callSet <- GenotypeCallXdrSet$fromFiles(path);
      verbose && cat(verbose, "Genotype-call set:");
      verbose && print(verbose, callSet);
      callList[[kk]] <- callSet;
      rm(callSet);
    }
  } # for (kk ...)

  plot(glad, arrays=arrays, chromosomes=chromosomes, skip=skip, callList=callList, imageFormat=imageFormat, verbose=less(verbose));
  verbose && exit(verbose);
}

if (5 %in% methods) {
  verbose && enter(verbose, "Writing GLAD regions to tabular file");
  # Filter out regions with only little change
  nbrOfSnps <- c(3,Inf);
  smoothing <- c(-Inf,-0.15,+0.15,+Inf);
  writeRegions(glad, arrays=arrays, chromosomes=chromosomes, smoothing=smoothing, nbrOfSnps=nbrOfSnps, oneFile=TRUE, skip=FALSE, verbose=less(verbose));
  verbose && exit(verbose);
}

if (6 %in% methods) {
  verbose && enter(verbose, "Writing WIG tracks");
  writeWig(glad, arrays=arrays, chromosomes=chromosomes, verbose=less(verbose));
  verbose && exit(verbose);
}

gc();
verbose && exit(verbose);

############################################################################
# HISTORY:
# 2006-11-xx
# o Recreated.
############################################################################
