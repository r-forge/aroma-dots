for (zzz in 0) {

verbose <- Arguments$getVerbose(TRUE)

# Look for the 100K Hind data directory
path <- filePath("cel/Hind/", expandLinks="any")
if (!isDirectory(path))
  break


name <- "ceph-Hind01-90"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Raw CEL dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("ceph")) {
  ceph <- AffymetrixDataset$fromFiles(path=path)
  if (nbrOfArrays(ceph) > 4)
    ceph <- extract(ceph, 1:4)
}
print(ceph)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Raw APD dataset (read optimized)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("apd")) {
  filename <- paste(name, "apd", class(ceph)[1], sep=".")

  # Load existing dataset?
  if (file.exists(filename)) {
    apd <- AffymetrixDataset$load(filename)
  } else {
    mkdirs("raw/")
    apd <- AffymetrixDataset$fromFiles(path="raw")
    # Create new data set?
    if (nbrOfArrays(apd) == 0) {
      # Optimized for reading units
      verbose && enter("Creating optimized reading and writing map")
      writeMap <- readCdfUnitsWriteMap(ceph[[1]])
      readMap <- invertMap(writeMap);
      verbose && exit(verbose)
      apd <- writeApd(ceph, path="raw", writeMap=writeMap, verbose=verbose)
      setMap(apd, readMap)
      rm(readMap, writeMap)
    } else if (nbrOfArrays(apd) != nbrOfArrays(ceph)) {
      throw(" number of arrays.")
    }
    save(apd, filename)
  }
  gc()
}
print(apd)


# Assert correctness
if (!exists("t2")) {
  probesets <- 55 + 1:1000
  rho <- (nbrOfProbesets(ceph)/length(probesets));

  verbose && enter(verbose, "Reading ", length(probesets), 
                                               " probesets from CEL files")
  t1 <- t <- system.time({
    ps1 <- getProbesets(ceph, probesets=probesets)
  })
  t <- t / nbrOfArrays(ceph);
  verbose && cat(verbose, 
              sprintf("Total time: %.2fs => %.2fs for all %d probesets", t[3], rho*t[3], nbrOfProbesets(ceph)))
  verbose && exit(verbose)


  verbose && enter(verbose, "Reading ", length(probesets), 
                                     " probesets from APD files")
  t2 <- t <- system.time({
    ps2 <- getProbesets(apd, probesets=probesets, readMap=NULL)
  })
  t <- t / nbrOfArrays(apd);
  verbose && cat(verbose, 
              sprintf("Total time: %.2fs => %.2fs for all %d probesets", t[3], rho*t[3], nbrOfProbesets(ceph)))
  verbose && exit(verbose)

  verbose && enter(verbose, "Reading ", length(probesets), 
                                     " probesets from optimized APD files")
  t3 <- t <- system.time({
    ps3 <- getProbesets(apd, probesets=probesets)
  })
  t <- t / nbrOfArrays(apd);
  verbose && cat(verbose, 
              sprintf("Total time: %.2fs => %.2fs for all %d probesets", t[3], rho*t[3], nbrOfProbesets(ceph)))
  verbose && exit(verbose)

  verbose && cat(verbose, 
                sprintf("CEL/APD time ratio: %.2f", t1[3]/t2[3]))
  verbose && cat(verbose, 
                sprintf("APD/optimized APD time ratio: %.2f", t2[3]/t3[3]))
  verbose && cat(verbose, 
                sprintf("CEL/optimized APD time ratio: %.2f", t1[3]/t3[3]))

  stopifnot(identical(ps1, ps3))

  # Clean up
  rm(probesets,t)
#  rm(ps1,ps2,ps3)
  gc()
}


} # for (zzz in 0)
