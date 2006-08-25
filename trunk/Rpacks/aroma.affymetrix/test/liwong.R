source("init.R")

cdfStrandJoiner <- function(units) {
  if (is.list(units[[1]])) {
    applyCdfGroups(units, cdfMergeStrands);
  } else {
    lapply(units, FUN=unique);
  }
}

verbose <- Arguments$getVerbose(TRUE);


ds <- AffymetrixCelSet$fromFiles("cel/Hind/");
print(ds);

cdf <- getCdf(ds);
setRestructor(cdf, cdfStrandJoiner);

model <- AffymetrixLiWongModel(ds);
print(model);

paf <- getProbeAffinities(model);
print(paf);

f <- fit(model, units=50000);
