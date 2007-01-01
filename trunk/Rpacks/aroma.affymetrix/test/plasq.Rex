library(aroma.affymetrix);
source("init.R");
source("ui.R");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
paths <- "probeData";
pattern <- "_(199|200)[0-9]";
if (exists("dataSets", mode="list")) {
  ds <- dataSets[[1]];
  name <- paste(getName(ds), paste(getTags(ds), collapse=","), sep=",");
  ans <- menu(c("<new data sets>", name), title="Select data set(s) to process:");
  if (ans == 1)
    dataSets <- selectDataSets(paths, pattern=pattern);
} else {
  dataSets <- selectDataSets(paths, pattern=pattern);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select chromosomes to fit
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# choices <- c(1:22,"X");
# if (!exists("chromosomes"))
#   chromosomes <- NULL;
# chromosomes <- selectMenu(choices=choices, selected=chromosomes, title="Select the chromosomes you want to fit (leave all blank to use all)");
# if (length(chromosomes) == 0) {
#   chromosomes <- choices;
# }

for (kk in seq(along=dataSets)) {
  ds <- dataSets[[kk]];
  ds <- extract(ds, 3:15);
  cdf <- getCdf(ds);

  plasq <- PlasqModel(ds);
}
