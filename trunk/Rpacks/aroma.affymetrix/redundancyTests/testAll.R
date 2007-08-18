files <- c(
# "test20070726,10k,CN.R",
# "test20070810,10k,profiling.R",
# "u133.test.mark.R",
# "test20070412,100K.R",
# "test20070412,100K,QN.R",
# "test20070412,100K,CN.R",
# "test20070625,100K,customCDF.R",
# "test20070809,100K,spatial.R",
 "test20070816,6.0.R"
)

for (file in files) {
  source(file, echo=TRUE);
  rm(list=setdiff(ls(), c("file", "files")));
  gc();
}
