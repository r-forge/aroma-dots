library(aroma.affymetrix);
#source("init.R");

source("refsets.init.R");

# Print distribution of dates
# dates <- lapply(dsR, FUN=function(ds) {
#   dates <- format(getTimestamps(ds), "%Y-%m-%d");
#   table(dates);
# })
# print(dates);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Model test set + complete reference set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ds <- clone(dsT);
append(ds, dsRall);
setTags(ds, c(getTags(ds), getName(dsRall)));
tags <- paste(getTags(ds), collapse=",");

testArrays <- which(getNames(ds) %in% getNames(dsT));
refArrays <- setdiff(seq(ds), testArrays);

rm(dsT, dsRall);
gc();

# Chromosomes for which GLAD should be fitted
allChromosomes <- c(rev(1:22), "X");
allChromosomes <- c(1:22, "X");
allChromosomes <- c(8,9,10,11,13,14);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup parallel jobs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
nbrOfRefArrays <- ceiling(seq(from=0, to=sqrt(length(refArrays)), length.out=12)^2);
#for (ff in nbrOfRefArrays) cat(file=sprintf("jobs/todo/%02d", ff), ff);
print(nbrOfRefArrays);

mkdirs("jobs/todo/");
mkdirs("jobs/running/");
mkdirs("jobs/done/");
mkdirs("jobs/progress/");


# Check remaining jobs
while (TRUE) {
  todo <- sort(list.files(path="jobs/todo", full.names=TRUE));
  if (length(todo) == 0) {
    verbose && cat(verbose, "No more jobs available. DONE!");
    break;
  }

  # Grab a random sample
  if (length(todo) > 1) {
    job <- sample(todo, size=1);
    job <- todo[1];
  } else {
    job <- todo;
  }
  name <- basename(job);

  # Try to obtain it
  verbose && enter(verbose, "Trying to obtain job ", name);
  runJob <- filePath("jobs/running", name);
  res <- file.rename(job, runJob);
  # If failed, try another one.
  if (!res) {
    verbose && exit(verbose, suffix="failed. Will try another one soon...");
    Sys.sleep(4*runif(n=1));
    next;
  }
#  cat(filerunJob, append=TRUE), " HOST=", System$getHostname());
  pb <- FileProgressBar(sprintf("jobs/progress/%s.progress", name));
  verbose && exit(verbose);

  tryCatch({
    # Grab the reference arrays in reverse order (first the most recent ones and so on).
    nbrOfRefs <- as.integer(name);
    arrays <- c(testArrays, rev(refArrays)[seq(length=nbrOfRefs)]);
    dsSub <- extract(ds, arrays);
    setTags(dsSub, c(getTags(dsSub), "time", sprintf("nref%02d", nbrOfRefs)));
    print(dsSub);
  
    # Quantile normalization, fitting PLM, fragment-length normalization
    cesFL <- estimateTotalCn(dsSub, progress=pb, verbose=verbose); # 37 steps
    rm(dsSub);
  
    # Define the GLAD model
    increase(pb);
    glad <- GladModel(cesFL);
    rm(cesFL);
    verbose && print(verbose, glad);
    increase(pb);
    
    # Get genome information
    cdf <- getCdf(glad);
    gi <- getGenomeInformation(cdf);
    rm(cdf);
    increase(pb);
  
    pbPerChr <- (pb$max - pb$value - 1) / length(allChromosomes);
    for (chr in allChromosomes) {
      verbose && enter(verbose, "Chromosome ", chr);
    
      units <- getUnitIndices(gi, chromosome=chr);
      nunits <- length(units);
      increase(pb, floor(0.05*pbPerChr));
 
      verbose && enter(verbose, "Fitting GLAD");
      fit(glad, arrays=testArrays, chromosomes=chr, .retResults=FALSE, verbose=verbose);
      verbose && exit(verbose);
      increase(pb, floor(0.70*pbPerChr));
    
      verbose && enter(verbose, "Plotting GLAD");
      plot(glad, arrays=testArrays, chromosomes=chr, verbose=verbose);
      verbose && exit(verbose);
      increase(pb, floor(0.25*pbPerChr));
    
#      verbose && enter(verbose, "Writing GLAD regions to tabular file");
#      # Filter out regions with only little change
#      smoothing <- c(-Inf,-0.15,+0.15,+Inf);
#      writeRegions(glad, arrays=testArrays, chromosomes=chr, smoothing=smoothing, 
#                                                            oneFile=TRUE, verbose=verbose);
#      verbose && exit(verbose);
    
      verbose && exit(verbose);
    }
  }, error=function(ex) {
    print(ex);
  })
  rm(gi,glad);
  setValue(pb, pb$max);
  rm(pb);  

  verbose && enter(verbose, "Trying to move job (", name, ") to done");
  doneJob <- filePath("jobs/done", name);
  res <- file.rename(runJob, doneJob);
  if (!res) {
    verbose && exit(verbose, suffix="failed.");
    Sys.sleep(4*runif(n=1));
    next;
  }
  verbose && exit(verbose);

  gc();
} # while(TRUE)

