
library("R.utils")

for (zzz in 0) {

verbose <- Arguments$getVerbose(TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. Scan for existing CEL files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# a) Scan for CEL files
path <- filePath("cel/", expandLinks="any")
files <- list.files(path=path, pattern="[.](c|C)(e|E)(l|L)$", full.names=TRUE)
files <- files[!file.info(files)$isdir]
if (length(files) == 0)
  break

cat("Converting CEL file to APD file: ", files[1], "\n")
apdfile <- celToApd(files[1])
cat("Created APD file: ", apdfile, "\n")
file.remove(apdfile)

cat("Create an optimal read map for CEL file: ", files[1], "\n")
cdffile <- findCdf(readCelHeader(files[1])$chiptype)
res <- cdfToApdMap(cdffile)

cat("Converting CEL file to APD file with an optimized read map: ", files[1], "\n")
apdfile <- celToApd(files[1], mapType=res$mapType)
cat("Created APD file: ", apdfile, "\n")

writeMap <- invertMap(res$readMap)
for (file in files[-1]) {
  cat("Converting CEL file to APD file with an optimized read map: ", file, "\n")
  apdfile <- celToApd(file, mapType=res$mapType, writeMap=writeMap)
  cat("Created APD file: ", apdfile, "\n")
}


} # for (zzz in 0)

rm(zzz, verbose, pattern)