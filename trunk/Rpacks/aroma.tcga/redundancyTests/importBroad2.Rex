if (interactive()) savehistory();
library("aroma.tcga");
sourceDirectory("R/", modifiedOnly=TRUE);
rootPath <- "tcgaData";

dataSet <- c("broad", "hudsonalpha", "harvard", "mskcc,A", "mskcc,B")[1];
paths <- getTcgaPaths(dataSet);

for (pp in seq(along=paths)) {
  path <- file.path(rootPath, paths[pp]);
  ds <- NULL;
  tryCatch({
    ds <- GenericDataFileSet$fromFiles(path, pattern=".*[.]data[.]txt$");
  }, error = function(ex) {}
  );
  if (is.null(ds)) {
    next;
  }
  print(ds);

  # Infer parameters
  archive <- basename(getPath(ds));
  dataType <- basename(getParent(getPath(ds), 1));
  dataType <- match.arg(dataType, c("cna", "snp"));
  platform <- basename(getParent(getPath(ds), 2));
  center <- basename(getParent(getPath(ds), 3));
  tumorType <- basename(getParent(getPath(ds), 5));
  
  # Identify unit names file
  if (platform == "genome_wide_snp_6") {
    require("aroma.affymetrix") || throw("Package not loaded: aroma.affymetrix");
    unf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6");
    dataTypes <- c("total", "fracB");
  } else if (platform == "") {
    dataTypes <- c("total", "fracB");
  } else if (platform == "") {
    dataTypes <- c("log2ratio,total");
  } else if (platform == "") {
    dataTypes <- c("log2ratio,total");
  } else {
    throw("Unknown TCGA platform: ", platform);
  }
  print(unf);
  chipType <- getChipType(unf, fullname=FALSE);
  
  rootPath <- "totalAndFracBData";
  for (kk in seq(ds)) {
    df <- getFile(ds, kk);
    db <- MageTabDataMatrixFile(getPathname(df));
    print(db);
    dataSet <- paste(c(archive, tumorType), collapse=",");
    path <- file.path(rootPath, dataSet, chipType);
  
    units <- NULL;
    for (cc in seq(from=2, to=nbrOfColumns(db))) {
      columnName <- getColumnNames(db)[cc];
      parts <- strsplit(columnName, split=",", fixed=TRUE)[[1]];
      parts <- setdiff(parts, "Signal");
      fullname <- paste(parts, collapse=",");
    
      fullnames <- paste(fullname, dataTypes, sep=",");
      filenames <- sprintf("%s.asb", fullnames);
      pathnames <- sapply(filenames, FUN=function(filename) {
        Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
      });
    
      # Nothing to do?
      isFiles <- sapply(pathnames, FUN=isFile);
      if (all(isFiles)) {
        print(pathnames);
        next;
      }
     
      print(pathnames);
    
      if (is.null(units)) {
        # Read unit (and allele) names
        unitNames <- readColumns(db, column=1)[,1];
        str(unitNames);
        if (dataType == "snp") {
          alleles <- rep(as.character(NA), length(unitNames));
          len <- nchar(unitNames);
          tail <- substr(unitNames, start=len-1, stop=len);
          for (aa in c("A", "B")) {
            idxs <- whichVector(tail == sprintf("-%s", aa));
            alleles[idxs] <- aa;
            unitNames[idxs] <- substr(unitNames[idxs], start=1, stop=len[idxs]-2);
            rm(idxs);
          }
        }
        str(unitNames);
        
        # Map unit names to unit indices
        units <- indexOf(unf, names=unitNames);
        str(units);
        print(summary(units));
        
        # Sanity check
        if (anyMissing(units)) {
          missing <- unitNames[is.na(units)];
          throw("There are ", length(missing), " unknown unit names: ", 
                                   paste(head(missing, 3), collapse=", "));
        }
      
        rm(unitNames);
      } # if (is.null(units))
    
      # Read input data
      values <- readColumns(db, column=cc, colClass="double")[,1];
      str(values);
    
      srcFile <- list(
        filename = getFilename(db),
        checksum = getChecksum(db),
        column = cc,
        columnName = columnName,
        valuesChecksum = digest(values)
      );
    
      # Transform to total or (total,beta)
      naValue <- as.double(NA);
      if (dataType == "snp") {
        data <- matrix(naValue, nrow=nbrOfUnits(unf), ncol=2);
        idxs <- whichVector(is.na(alleles));
        data[units[idxs],1] <- values[idxs];
        idxs <- whichVector(alleles == "A");
        data[units[idxs],1] <- values[idxs];
        idxs <- whichVector(alleles == "B");
        data[units[idxs],1] <- data[units[idxs],1] + values[idxs];
        data[units[idxs],2] <- values[idxs] / data[units[idxs],1];
        rm(idxs);
      } else if (dataType == "cna") {
        data <- matrix(naValue, nrow=nbrOfUnits(unf), ncol=1);
        data[units,1] <- values;
      }
      rm(values, units);

      # Write total
      for (dd in whichVector(!isFiles)) {
        pathname <- pathnames[dd];
        df <- AromaUnitTotalCnBinaryFile$allocateFromUnitNamesFile(unf, filename=pathname, path=NULL);
        footer <- readFooter(df);
        footer$srcFile <- srcFile;
        writeFooter(df, footer);
        df[,1] <- data[,dd];
        print(df);
        rm(df);
      } # for (dd ...)
      rm(data);

    } # for (cc ...)

  } # for (kk ...)

} # for (pp ...)
