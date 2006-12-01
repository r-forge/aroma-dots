library(aroma.affymetrix);
source("init.R");

source("refsets.init.R");

# Print distribution of dates
dates <- lapply(dsR, FUN=function(ds) {
  dates <- format(getTimestamps(ds), "%Y-%m-%d");
  table(dates);
})
print(dates);


ans <- trim(readline("Enter 1 or 2: "));
dsR0 <- dsR[as.integer(ans)];

for (kk in seq(along=dsR0)) {
  ds <- clone(dsT);
  append(ds, dsR0[[kk]]);
  setTags(ds, c(getTags(ds), getName(dsR0[[kk]])));
  tags <- paste(getTags(ds), collapse=",");

  # Quantile normalization, fitting PLM, fragment-length normalization
  cesFL <- estimateTotalCn(ds, verbose=verbose);

  # Define the GLAD model
  glad <- GladModel(cesFL);
  verbose && print(verbose, glad);

  # Get genome information
  cdf <- getCdf(glad);
  gi <- getGenomeInformation(cdf);

  allChromosomes <- c(rev(1:22), "X");
  allChromosomes <- c(1:22, "X");
  allChromosomes <- allChromosomes[1];
  for (chr in allChromosomes) {
    verbose && enter(verbose, "Chromosome ", chr);

    units <- getUnitIndices(gi, chromosome=chr);
    nunits <- length(units);

    verbose && enter(verbose, "Fitting GLAD");
    fit(glad, arrays=1:3, chromosomes=chromosome, .retResults=FALSE, verbose=verbose);
    verbose && exit(verbose);

    verbose && enter(verbose, "Plotting GLAD");
    plot(glad, arrays=1:3, chromosomes=chromosome, verbose=verbose);
    verbose && exit(verbose);

    verbose && enter(verbose, "Writing GLAD regions to tabular file");
    # Filter out regions with only little change
    smoothing <- c(-Inf,-0.15,+0.15,+Inf);
    writeRegions(glad, arrays=1:3, chromosomes=chromosome, smoothing=smoothing, 
                                                        oneFile=TRUE, verbose=verbose);
    verbose && exit(verbose);

    verbose && exit(verbose);
  }

  verbose && exit(verbose);
}
