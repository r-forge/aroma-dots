if (interactive()) savehistory();
library("aroma.tcga");
library("R.menu");
sourceDirectory("R/", modifiedOnly=TRUE);
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

rootPath <- "rawData/";
rootPath <- Arguments$getReadablePath(rootPath);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify all data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pattern <- TcgaDccData$getDataSetPatterns()$dataset;
dataSets <- list.files(path=rootPath, pattern=pattern);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select by center
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
centers <- gsub(pattern, "\\1", dataSets);
uCenters <- sort(unique(centers));
if (interactive()) {
  center <- textMenu(uCenters, value=TRUE);
} else {
  center <- centers[1];
}
verbose && cat(verbose, "TCGA center:");
verbose && print(verbose, center);
dataSets <- dataSets[centers == center];


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select chip type
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chipTypeList <- lapply(dataSets, FUN=function(dataSet) {
  path <- file.path(rootPath, dataSet);
  dirs <- list.files(path=path);
  dirs <- gsub("[.](lnk|LNK)", "", dirs);
  dirs;
});
uChipTypes <- sort(unique(unlist(chipTypeList, use.names=FALSE)));
if (interactive()) {
  chipType <- textMenu(uChipTypes, value=TRUE);
} else {
  chipType <- uChipTypes[1];
}
verbose && cat(verbose, "Chip type:");
verbose && print(verbose, chipType);
keep <- sapply(chipTypeList, FUN=function(chipTypes) { 
  any(chipTypes == chipType);
});
dataSets <- dataSets[keep];


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Select type of data to import
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (center == "broad.mit.edu") {
  typeList <- list(
    "ismpolish" = list(
      clazz = BroadAlleleSpecificThetaTcgaDataFileSet
    ),
    "5NN" = list(
      clazz = BroadTotalCopyNumberTcgaDataFileSet,
      exportFcn = exportTotalAndFracB
    ),
    "ASCN" = list(
      clazz = BroadAlleleSpecificCopyNumberTcgaDataFileSet
    )
  );
} else if (center == "hudsonalpha.org") {
  typeList <- list(
    "B_Allele_Freq" = list(
      clazz = HudsonAlphaBAlleleFreqTcgaDataFileSet,
      exportFcn = exportTotalAndFracB
    )
  );
} else {
  throw("Unsupported center: ", center);
}
if (interactive()) {
  typeName <- textMenu(names(typeList), value=TRUE);
} else {
  typeName <- names(typeList)[2];
}
verbose && cat(verbose, "Type of data: ", typeName);
type <- typeList[[typeName]];

clazz <- type$clazz;
exportFcn <- type$exportFcn;
if (is.null(exportFcn)) {
  throw("No export function specified/implemented: ", typeName);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get unit names file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Get unit names file");
if (chipType == "GenomeWideSNP_6") {
  require("aroma.affymetrix") || throw("Package not loaded: aroma.affymetrix");
  unf <- AffymetrixCdfFile$byChipType("GenomeWideSNP_6");
} else if (chipType == "Human1M-Duo") {
  unf <- TextUnitNamesFile$byChipType(chipType, tags="TCGA,OV,hudsonalpha.org,HB20090506");
} else if (chipType == "HG-CGH-244A") {
  unf <- TextUnitNamesFile$byChipType(chipType, tags="HB20090419");
} else if (chipType == "HG-CGH-415K") {
  unf <- TextUnitNamesFile$byChipType(chipType, tags="adf,TCGA,hms.harvard.edu,HB20090824");
} else {
  throw("Unknown TCGA chip type: ", chipType);
}
verbose && print(verbose, unf);
verbose && exit(verbose);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (ii in seq(along=dataSets)) {
  dataSet <- dataSets[ii];
  verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                            ii, dataSet, length(dataSets)));

  ds <- clazz$byName(dataSet, subdirs=chipType, paths=rootPath);
  setFullName(ds, dataSet);
  verbose && print(verbose, ds);
  # Nothing to do?
  if (nbrOfFiles(ds) == 0) {
    verbose && cat(verbose, "No TCGA data files. Nothing to import.");
    verbose && exit(verbose);
    next;  
  }

  # Import each data file
  verbose && enter(verbose, "Exporting TCGA data as aroma files");
  dsList <- exportFcn(ds, unf=unf, verbose=less(verbose, 5));
  verbose && print(verbose, dsList);
  verbose && exit(verbose);
  
  verbose && exit(verbose);
} # for (ii ...)

stop();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data specific classes
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

signalPattern <- ",Signal(|_[AB])$";

for (pp in seq(along=paths)) {
  rootPath <- "tcgaData";
  path <- file.path(rootPath, paths[pp]);

  verbose && enter(verbose, sprintf("Path #%d of %d", pp, length(paths)));
  verbose && cat(verbose, "Path: ", path);

  # Infer parameters
  archive <- basename(getPath(ds));
  dataType <- basename(getParent(getPath(ds), 1));
  dataType <- match.arg(dataType, c("cna", "snp"));
  platform <- basename(getParent(getPath(ds), 2));
  center <- basename(getParent(getPath(ds), 3));
  tumorType <- basename(getParent(getPath(ds), 5));

  verbose && cat(verbose, "Platform: ", platform);

  verbose && enter(verbose, "Identifying unit names file");
  # Identify unit names file
  if (platform == "genome_wide_snp_6") {
    if (type == "5NN") {
      dataTypes <- "total";
    } else {
      dataTypes <- c("total", "fracB");
    }
  } else if (platform == "human1mduo") {
    dataTypes <- c("total", "fracB");
  } else if (platform == "hg-cgh-244a") {
    dataTypes <- c("log2ratio,total");
  } else if (platform == "hg-cgh-415k_g4124a") {
    dataTypes <- c("log2ratio,total");
  } else {
    throw("Unknown TCGA platform: ", platform);
  }
  verbose && print(verbose, unf);
  verbose && exit(verbose);

  rootPath <- "totalAndFracBData";
  for (kk in seq(ds)) {
    df <- getFile(ds, kk);
    verbose && enter(verbose, sprintf("Data file #%d of %d", kk, length(ds)));
    verbose && cat(verbose, "Full name: ", getFullName(df));

    db <- clazz(getPathname(df));
    verbose && print(verbose, db);

    dataSet <- paste(c(archive, tumorType, tags), collapse=",");
    verbose && cat(verbose, "Output data set: ", dataSet);

    path <- file.path(rootPath, dataSet, chipType);
    verbose && cat(verbose, "Output path: ", path);
  
    verbose && enter(verbose, "Reading column names");
    columnNames <- getColumnNames(df, verbose=less(verbose,5));
    verbose && exit(verbose);

    units <- NULL;
    signalColumnIdxs <- grep(signalPattern, columnNames);
    for (aa in seq(along=signalColumnIdxs)) {
      cc <- signalColumnIdxs[aa];
      columnName <- columnNames[cc];
      verbose && enter(verbose, sprintf("Array #%d ('%s') of %d",
                                    aa, columnName, length(signalColumnIdxs)));

      parts <- strsplit(columnName, split=",", fixed=TRUE)[[1]];
      parts <- setdiff(parts, "Signal");
      parts <- setdiff(parts, "Log2(Rtum/Rnorm)");
      parts <- setdiff(parts, "signal:Log2");
      parts <- setdiff(parts, "normalizedLog2Ratio");
      parts <- setdiff(parts, "Log2Ratio");
      fullname <- paste(parts, collapse=",");
    
      fullnames <- paste(fullname, dataTypes, sep=",");
      filenames <- sprintf("%s.asb", fullnames);
      pathnames <- sapply(filenames, FUN=function(filename) {
        Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
      });
    
      # Nothing to do?
      isFiles <- sapply(pathnames, FUN=isFile);
      if (all(isFiles)) {
        verbose && cat(verbose, "All data is already imported");
        verbose && print(verbose, pathnames);
        verbose && exit(verbose);
        next;
      }
     
      verbose && print(verbose, pathnames);

      if (is.null(units)) {
        verbose && enter(verbose, "Identifying unit indices");

        verbose && enter(verbose, "Reading unit and allele names");
        # Read unit (and allele) names
        unitNames <- readColumns(db, column=1)[,1];
        verbose && str(verbose, unitNames);
        verbose && exit(verbose);

        if (dataType == "snp" && type == "ismpolish") {
          verbose && enter(verbose, "Identifying (unit, allele) names");
          alleles <- rep(as.character(NA), length(unitNames));
          len <- nchar(unitNames);
          tail <- substr(unitNames, start=len-1, stop=len);
          for (bb in c("A", "B")) {
            idxs <- whichVector(tail == sprintf("-%s", bb));
            alleles[idxs] <- bb;
            unitNames[idxs] <- substr(unitNames[idxs], start=1, stop=len[idxs]-2);
            rm(idxs);
          }
          verbose && cat(verbose, "Unit names:");
          verbose && str(verbose, unitNames);
          verbose && cat(verbose, "Alleles:");
          verbose && str(verbose, alleles);
          verbose && exit(verbose);
        }

        verbose && enter(verbose, "Mapping unit names to unit indices");
        # Map unit names to unit indices
        units <- indexOf(unf, names=unitNames);
        verbose && str(verbose, units);
        verbose && summary(verbose, units);
        verbose && exit(verbose);

        # Sanity check
        if (anyMissing(units)) {
          missing <- unitNames[is.na(units)];
          throw("There are ", length(missing), " unknown unit names: ", 
                                   paste(head(missing, 3), collapse=", "));
        }
      
#        rm(unitNames);
        verbose && exit(verbose);
      } # if (is.null(units))

      # Read input data
      verbose && enter(verbose, "Reading data");
      colClass <- "double";
      values <- readColumns(db, column=cc, colClass=colClass)[,1];
      verbose && str(verbose, values);
      verbose && summary(verbose, values);
      verbose && exit(verbose);

      # Transform to total or (total,beta)
      verbose && enter(verbose, "Restructuring and transforming signals on file");
      naValue <- as.double(NA);
      if (dataType == "snp" && chipType == "Human1M-Duo") {
        data <- matrix(naValue, nrow=nbrOfUnits(unf), ncol=2);
        data[units,1] <- (values[,1] + values[,2]);
        data[units,2] <- values[,2] / data[units,1];
      } else if (dataType == "snp" && chipType == "GenomeWideSNP_6" && type != "5NN") {
        verbose && enter(verbose, "Transforming to {(C,beta)}");
        data <- matrix(naValue, nrow=nbrOfUnits(unf), ncol=2);
        idxs <- whichVector(is.na(alleles));
        data[units[idxs],1] <- values[idxs];
        idxs <- whichVector(alleles == "A");
        data[units[idxs],1] <- values[idxs];
        idxs <- whichVector(alleles == "B");
        data[units[idxs],1] <- data[units[idxs],1] + values[idxs];
        data[units[idxs],2] <- values[idxs] / data[units[idxs],1];
        rm(idxs);
        verbose && exit(verbose);
      } else if (dataType == "cna" || type == "5NN") {
        verbose && enter(verbose, "Transforming to {C}");
        data <- matrix(naValue, nrow=nbrOfUnits(unf), ncol=1);
        data[units,1] <- values;
        verbose && exit(verbose);
      } else {
        throw("Unknown data format.");
      }
#      rm(values);
      verbose && str(verbose, data);
      verbose && summary(verbose, data);
      verbose && exit(verbose);

      # Write total
      for (dd in whichVector(!isFiles)) {
        verbose && enter(verbose, sprintf("Writing output file #%d of %d", dd, sum(!isFiles)));

        pathname <- pathnames[dd];
        verbose && cat(verbose, "Output pathname: ", pathname);

        verbose && enter(verbose, "Allocating empty file");
        df <- AromaUnitTotalCnBinaryFile$allocateFromUnitNamesFile(unf, 
                                          filename=pathname, path=NULL);
        footer <- readFooter(df);
        footer$srcFile <- srcFile;
        writeFooter(df, footer);
        verbose && exit(verbose);

        verbose && enter(verbose, "Write signals");
        df[,1] <- data[,dd];
        verbose && exit(verbose);

        verbose && print(verbose, df);
#        rm(df);
        verbose && exit(verbose);
      } # for (dd ...)
#      rm(data);

      verbose && exit(verbose);
    } # for (aa ...)

    verbose && exit(verbose);
  } # for (kk ...)

  verbose && exit(verbose);
} # for (pp ...)


############################################################################
# HISTORY:
# 2009-09-24
# o Added verbose output.
# 2009-xx-xx
# o Created.
############################################################################
