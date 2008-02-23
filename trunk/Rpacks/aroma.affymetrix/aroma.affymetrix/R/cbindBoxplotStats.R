setMethodS3("cbindBoxplotStats", "list", function(stats, ...) {
  # Append 'group' stats
  for (kk in seq(along=stats)) {
    stats[[kk]][["group"]] <- rep(kk, times=length(stats[[kk]][["out"]]));
  }

  # Restructure so it is recognized by graphics::bxp().
  bxpStats <- list();

  for (field in names(stats[[1]])) {
    args <- lapply(stats, FUN=.subset2, field);
    value <- do.call("c", args);
    if (field == "stats") {
      value <- matrix(value, nrow=5);
    } else if (field == "conf") {
      value <- matrix(value, nrow=2);
    }
    bxpStats[[field]] <- value;
  }

  bxpStats[["names"]] <- names(stats);

  bxpStats;
}, protected=TRUE)


##########################################################################
# HISTORY:
# 2008-02-22
# o Created.
##########################################################################

