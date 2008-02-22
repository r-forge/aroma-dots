setMethodS3("getUnitGroupNamesFromUgcMap", "AffymetrixCdfFile", function(this, ugcMap, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting unit and group names from CDF");

  units <- ugcMap[,"unit"];

  # Get unit names
  res <- data.frame(
    unitName = getUnitNames(this, units=units),
    stringsAsFactors = FALSE
  );
  verbose && str(verbose, res);

  verbose && enter(verbose, "Reading all group names for units of interest");
  uniqueUnits <- unique(units);
  groupNames <- readCdfGroupNames(getPathname(this), units=uniqueUnits);
  verbose && cat(verbose, "First unit:");
  verbose && str(verbose, groupNames[1]);
  verbose && exit(verbose);

  groups <- ugcMap[,"group"];
  verbose && enter(verbose, "Build (unit name, group name) map");
  verbose && cat(verbose, "Number of units: ", length(uniqueUnits));
  for (kk in seq(along=uniqueUnits)) {
    if (kk %% 100 == 0) {
      verbose && writeRaw(verbose, length(uniqueUnits)-kk, ", ");
    }

    # Find matching rows
    ### rr <- which(units %in% .uniqueUnits[kk]);
    rr <- which(units %in% .subset(uniqueUnits, kk));
    # All groups for this unit
    ### gg <- groups[rr];
    ### res[rr,"groupName"] <- groupNames[[kk]][gg];
    gg <- .subset(groups, rr);
    res[rr,"groupName"] <- .subset(.subset2(groupNames, kk), gg);
  }
  verbose && cat(verbose, "0.");
  verbose && exit(verbose);

  if (nrow(res) != nrow(ugcMap)) {
    throw("Internal error: Number of extract unit and group names does not match the number of rows in the UGC map: ", nrow(res), " != ", nrow(ugcMap));
  }

  verbose && exit(verbose);

  res;
}) # getUnitGroupNamesFromUgcMap()


############################################################################
# HISTORY:
# 2008-02-05
# o Created.
############################################################################
