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
  }

  # Remove duplicates
  ds <- extract(ds, !isDuplicated(ds));

  # Order by timestamps
  o <- order(as.POSIXct(getTimestamps(ds)));
  ds <- extract(ds, o);

  setName(ds, "AGRFall");

  dsRall <- ds;
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Grouped reference sets (by date)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("dsR", mode="list")) {
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
    ds0 <- extract(ds, arrays)
    setName(ds0, sprintf("AGRF%02d", gg));
    dsR <- c(dsR, list(ds0));
  }
}


setMethodS3("estimateTotalCn", "AffymetrixCelSet", function(ds, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
  verbose && print(verbose, normQ);
  dsN <- process(normQ, verbose=verbose);
  verbose && exit(verbose);

  # Fit PLM
  verbose && enter(verbose, "Summarazing probe signals");
  plm <- RmaCnPlm(dsN, mergeStrands=TRUE, combineAlleles=TRUE);
  verbose && print(verbose, plm);
  fit(plm, moreUnits=3, verbose=verbose);
  verbose && exit(verbose);

  # Normalize for PCR fragment-length effects
  verbose && enter(verbose, "Normalizing summarized SNP signals");
  ces <- getChipEffects(plm);
  verbose && print(verbose, ces);
  normFL <- FragmentLengthNormalization(ces, subsetToFit=1/3);
  verbose && print(verbose, normFL);
  cesFL <- process(normFL, verbose=verbose);
  verbose && print(verbose, cesFL);
  verbose && exit(verbose);

  verbose && exit(verbose);

  cesFL;
}) # estimateTotalCn()
