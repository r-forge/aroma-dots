library(aroma.affymetrix);
source("init.R");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The test set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "raw/CarmichaelC_etal_2006-500k/Mapping250K_Nsp";
dsT <- AffymetrixCelSet$fromFiles(path);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The complete reference set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("dsRall")) {
  ds <- NULL;
  names <- c("SlaterH_etal_2006", "SinclairA_etal_2006", "AGRF_2006a-Anonymous");
  for (name in names) {
    path <- filePath("raw", name, "Mapping250K_Nsp");
    ds0 <- AffymetrixCelSet$fromFiles(path);
    if (is.null(ds)) {
      ds <- ds0;
    } else {
      append(ds, ds0);
    }
    rm(ds0);
  }
  gc();

  # Remove duplicates
  ds <- extract(ds, !isDuplicated(ds));

  # Order by timestamps
  o <- order(as.POSIXct(getTimestamps(ds)));
  ds <- extract(ds, o);

  setName(ds, "AGRFall");

  dsRall <- ds;
  rm(ds); gc();
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Grouped reference sets (by date)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (FALSE && !exists("dsR", mode="list")) {
  ds <- dsRall;

  # Split into subgroups
  dates <- format(getTimestamps(ds), "%Y-%m-%d");
  print(table(dates));
  # 2006-05-16 2006-07-13 2006-07-14 2006-09-29 2006-10-04 2006-11-01
  #          9          5          5          6          5          5
  
  dsR <- list();
  patterns <- c("2006-05-16", "2006-07-(13|14)", "2006-(10-04|09-29)");
  for (gg in seq(along=patterns)) {
    pattern <- patterns[gg];
    arrays <- grep(pattern, dates);
    arrays <- arrays[1:9];
    dates[arrays] <- NA;  # Avoid mistakes
    ds0 <- extract(ds, arrays);
    setName(ds0, sprintf("AGRF%02d", gg));
    dsR <- c(dsR, list(ds0));
    rm(ds0);
  }
  rm(ds);
  gc();
}


setMethodS3("estimateTotalCn", "AffymetrixCelSet", function(ds, ..., progress, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'progress':
  pb <- progress;

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Estimating SNP signals for total copy-number analysis");
  verbose && print(verbose, ds);

  # Normalize
  verbose && enter(verbose, "Normalizing probe signals");
  normQ <- QuantileNormalizer(ds, subsetToAvg=1/3);
  rm(ds);
  verbose && print(verbose, normQ);
  increase(pb, 1);
  dsN <- process(normQ, verbose=verbose);
  rm(normQ); 
  increase(pb, 10);
  gc();
  verbose && exit(verbose);

  # Fit PLM
  verbose && enter(verbose, "Summarazing probe signals");
  plm <- RmaCnPlm(dsN, mergeStrands=TRUE, combineAlleles=TRUE);
  increase(pb, 1);
  rm(dsN);
  verbose && print(verbose, plm);
  fit(plm, moreUnits=1, verbose=verbose);
  increase(pb, 20);
  ces <- getChipEffects(plm);
  rm(plm);
  gc();
  verbose && exit(verbose);

  # Normalize for PCR fragment-length effects
  verbose && enter(verbose, "Normalizing summarized SNP signals");
  verbose && print(verbose, ces);
  normFL <- FragmentLengthNormalization(ces, subsetToFit=1/3);
  increase(pb, 1);
  rm(ces);
  verbose && print(verbose, normFL);
  cesFL <- process(normFL, verbose=verbose);
  increase(pb, 4);
  rm(normFL);
  gc();
  verbose && print(verbose, cesFL);
  verbose && exit(verbose);

  verbose && exit(verbose);

  cesFL;
}) # estimateTotalCn()
