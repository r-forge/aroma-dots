###########################################################################/**
# @set "class=SnpChipEffectSet"
# @RdocMethod exportToAPT
#
# @title "Export probeset summaries to an Affymetrix Power Tool (APT) result file"
#
# \description{
#  @get "title" named *.summary.txt.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the APT result file.}
#   \item{path}{An optional path to the file.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns (invisible) the pathname of the generated file.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("exportToAPT", "SnpChipEffectSet", function(this, filename=sprintf("%s,%s.summary.txt", getFullName(this), getChipType(getCdf(this), fullname=FALSE)), path=NULL, sampleNames=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (is.null(sampleNames)) {
    sampleNames <- sapply(this, getFileName);
  }

  if (!this$mergeStrands) {
    throw("Cannot export chip effects to an APT summary file, because the strands are not merged, which is required by the APT applications.");
  }

  verbose && enter(verbose, "Exporting unit group summaries to an APT summary file");
  verbose && cat(verbose, "Pathname: ", pathname);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identifying units to be exported (those with two groups)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying units to be exported");  
  cdf <- getCdf(this);
  unitSizes <- getUnitSizes(cdf);
  units <- which(unitSizes >= 2);
  rm(unitSizes);
  verbose && cat(verbose, "Number of units to be exported: ", length(units));
  verbose && str(verbose, units);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying (unit, group, cell) map");  
  ugcMap <- getUnitGroupCellMap(cdf, units=units, verbose=less(verbose));
  rm(units);
  keep <- (ugcMap[,"group"] <= 2)
  ugcMap <- ugcMap[keep,,drop=FALSE];
  rm(keep);
  verbose && exit(verbose);

  # Get unit names
  unitGroupNames <- getUnitNames(cdf, units=ugcMap[,"unit"]);
  unitGroupNames <- paste(unitGroupNames, c("A", "B"), sep="-");

  # Get cells
  cells <- ugcMap[,"cell"];
  rm(ugcMap);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Open file & assert file format
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  con <- file(pathname, open="w");
  on.exit(close(con));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  guid <- "0000065391-1175407180-0000005910-0000025527-0000027480";
  execGuid <- "0000065535-1175407155-0000018162-0000000009-0000007239";
  execVersion <- "$Id: gtype-probeset-genotype.cpp,v 1.1.2.1 2006/11/21 05:40:33 sbanvi Exp $-1.5.6";
  cmd <- "gtype-probeset-genotype --verbose 2 --analysis quant-norm.sketch=50000.usepm=true,pm-only,sumz --out-dir data/out/nsp --cdf-file lib/cdf/Mapping250K_Nsp.cdf --no-chp-output --chrX-force --summaries --cel-files Affymetrix_2006-HapMap270.CEU.founders,files.txt";

  cdf <- getCdf(this);
  chipType <- gsub(",monocell", "", getChipType(cdf));
  header <- c(
    "guid" = guid,
    "exec_guid" = execGuid,
    "exec_version" = execVersion,
    "create_date" = as.character(Sys.time()),
    "cmd" = cmd,
    "chip_type" = chipType,
    "lib_set_name" = getPathname(cdf),
    "lib_set_version" = getPathname(cdf)
  );

  header <- paste(names(header), "=", header, sep="");
  header <- paste("#%", header, sep="");

  writeLines(con=con, header, sep="\n");
  cat(file=con, "############################################################\n");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Column names
  colnames <- c("probeset_id", sampleNames);
  colnames <- paste(colnames, collapse="\t");
  cat(file=con, colnames, "\n", sep="");

  bytesPerRow <- 64+2+1 + 8*nbrOfArrays(this);
  bytesPerChunk <- 25e6;
  cellsPerChunk <- floor(bytesPerChunk / bytesPerRow);
  nbrOfChunks <- ceiling(length(cells) / cellsPerChunk);

  count <- 1;
  while (length(cells) > 0) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", count, nbrOfChunks));
    # Get cell indices to be read
    if (length(cells) > cellsPerChunk) {
      head <- 1:cellsPerChunk;
    } else {
      head <- 1:length(cells);
    }

    data <- getData(this, indices=cells[head], fields="intensities")$intensities;
    data <- round(data, digits=5);  # As APT does!
    rownames(data) <- unitGroupNames[head];
    verbose && str(verbose, data);

    verbose && enter(verbose, "Appending data table to output file");
    write.table(file=con, data, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE);
    verbose && exit(verbose);

    rm(data);

    cells <- cells[-head];
    unitGroupNames <- unitGroupNames[-head];

    # Garbage collect
    gc <- gc();
    verbose && print(verbose, gc);

    count <- count + 1;

    verbose && exit(verbose);
  } # while(...)


  verbose && exit(verbose);

  invisible(pathname);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2007-04-01
# o Created.
############################################################################
