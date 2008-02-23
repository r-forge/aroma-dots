setMethodS3("extractMatrix", "ParameterCelFile", function(this, units=NULL, ..., field=c("intensities", "stdvs", "pixels"), returnUgcMap=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  cdf <- getCdf(this);

  # Argument 'units':
  if (is.null(units)) {
  } else {
    units <- Arguments$getIndices(units, range=c(1, nbrOfUnits(cdf)));
  }

  # Argument 'field':
  if (length(field) > 1)
    field <- field[1];

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting data for the array set");

  verbose && enter(verbose, "Getting unit-to-cell map");
  ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
  ugcMap <- subset(ugcMap, ...);
  verbose && exit(verbose);

  if (nrow(ugcMap) == 0)
    throw("Nothing to return.");

  df <- matrix(NA, nrow=nrow(ugcMap), ncol=1);
  colnames(df) <- getName(this);

  gc <- gc();
  verbose && print(verbose, gc);

  verbose && enter(verbose, "Retrieving array data");
  df[, 1] <- getDataFlat(this, units=ugcMap, fields=field, 
                                           verbose=less(verbose))[, field];
  verbose && exit(verbose);

  verbose && exit(verbose);

  if (returnUgcMap)
    attr(df, "unitGroupCellMap") <- ugcMap;

  df;
})


setMethodS3("extractDataFrame", "ParameterCelFile", function(this, addNames=FALSE, addUgcMap=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting data for the array");
  data <- extractMatrix(this, ..., returnUgcMap=TRUE, 
                                                   verbose=less(verbose, 1));

  ugcMap <- attr(data, "unitGroupCellMap");
  attr(data, "unitGroupCellMap") <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  if (addUgcMap) {
    verbose && enter(verbose, "Merging UGC map and extracted data");
    ugcMap <- as.data.frame(ugcMap);
    data <- cbind(ugcMap, data);

    if (addNames) {
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
    }

    verbose && exit(verbose);
  }

  if (addNames) {
    verbose && enter(verbose, "Appending unit and group names from CDF");
    cdf <- getCdf(this);
    verbose && cat(verbose, "CDF chip type: ", 
                                        getChipType(cdf, fullname=TRUE));
    ugNames <- getUnitGroupNamesFromUgcMap(cdf, ugcMap=ugcMap, 
                                              verbose=less(verbose, 10));
    rm(cdf, ugcMap);
    verbose && cat(verbose, "(unit, group) names: ");
    verbose && str(verbose, ugNames);

    ugNames <- as.data.frame(ugNames);
    data <- cbind(ugNames, data);
    rm(ugNames);

    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  data;
}) # extractDataFrame()


############################################################################
# HISTORY:
# 2007-02-22
# o Added extractDataFrame() for ParameterCelFile as well.
# o Generalized to ParameterCelFile and moved into aroma.affymetrix.
# 2008-02-11 [EP]
# o Created for ChipEffectFile and FirmaFile.
############################################################################ 
