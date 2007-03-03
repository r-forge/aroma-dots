setConstructorS3("AromaUgpFile", function(...) {
  extend(AromaGenomePositionFile(...), "AromaUgpFile")
})

setMethodS3("getFilenameExtension", "AromaGenomePositionFile", function(static, ...) {
  "ugp";
}, static=TRUE)

setMethodS3("nbrOfUnits", "AromaUgpFile", function(this, ...) {
  nbrOfElements(this, ...);
})


setMethodS3("getCdf", "AromaUgpFile", function(this, ...) {
  cdf <- this$.cdf;
  if (is.null(cdf)) {
    chipType <- getChipType(this);
    cdf <- AffymetrixCdfFile$fromChipType(chipType);
    this$.cdf <- cdf;
  }
  
  cdf;
})

setMethodS3("createFromCdf", "AromaUgpFile", function(static, cdf, ...) {
  chipType <- getChipType(cdf);
  create(static, chipType=chipType, nbrOfElements=nbrOfUnits(cdf), ...);
}, static=TRUE)


setMethodS3("indexOfElements", "AromaUgpFile", function(this, names, ...) {
  # Look up unit names from CDF
  cdf <- getCdf(this);
  idxs <- match(names, getUnitNames(cdf));
  idxs;
}, protected=TRUE)

setMethodS3("getUnitsAt", "AromaUgpFile", function(this, ...) {
  getElementsAt(this, ...);
})


setMethodS3("importFromGenomeInformation", "AromaUgpFile", function(this, gi, ..., verbose=FALSE) {
  if (!inherits(gi, "GenomeInformation")) {
    throw("Argument 'gi' is not a GenomeInformation object: ", class(gi)[1]);
  }

  # AD HOC patch, since units==NULL does not work./HB 2007-03-03
  units <- seq_len(nbrOfUnits(gi));
  data <- getData(gi, units=units, fields=c("chromosome", "physicalPosition"));

  chr <- data[,"chromosome"];
  if (is.character(chr)) {
    chr[chr == "X"] <- 23;
    chr[chr == "Y"] <- 23;
    suppressWarnings({
      chr <- as.integer(chr);
    })
  }
  
  pos <- data[,"physicalPosition"];
  suppressWarnings({
    pos <- as.integer(pos);
  })

  updateData(this, chromosome=chr, physicalPosition=pos);
})


setMethodS3("createFromGenomeInformation", "AromaUgpFile", function(static, gi, ..., verbose=FALSE) {
  if (!inherits(gi, "GenomeInformation")) {
    throw("Argument 'gi' is not a GenomeInformation object: ", class(gi)[1]);
  }

  chipType <- getChipType(gi);
  cdf <- AffymetrixCdfFile$fromChipType(chipType);
  ugp <- createFromCdf(cdf, ...);
  importFromGenomeInformation(ugp, gi);
}, static=TRUE)



############################################################################
# HISTORY:
# 2007-03-03
# o Now inherits from generic AromaGenomePositionFile.
# 2007-03-02
# o Created. Can import genome information data.
############################################################################
