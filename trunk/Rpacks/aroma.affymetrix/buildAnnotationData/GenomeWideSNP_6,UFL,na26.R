if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "26";
user <- "HB";
datestamp <- "20080722";

chipType <- "GenomeWideSNP_6";
cdfTags <- "Full";
nbrOfEnzymes <- 2;

footer <- list(
  createdOn = format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
  createdBy = list(
    fullname = "Henrik Bengtsson", 
    email = "hb@stat.berkeley.edu"
  ),
  srcFiles = list()
);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup required annotation files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("cdf")) {
  cdf <- AffymetrixCdfFile$byChipType(chipType, tags=cdfTags);
  rm(csvList);
}
print(cdf);

if (!exists("csvList", mode="list")) {
  csvList <- list();

  tagsList <- c(
      main=sprintf(   ".na%s", naVersion),
      cn=sprintf(".cn.na%s", naVersion)
  );

  for (key in names(tagsList)) {
    tags <- tagsList[[key]];
    pathname <- AffymetrixNetAffxCsvFile$findByChipType(chipType, tags=tags);
    if (isFile(pathname)) {
      csv <- AffymetrixNetAffxCsvFile(pathname);
      srcFile <- list(
        filename=getFilename(csv), 
        filesize=getFileSize(csv), 
        checksum=getChecksum(csv)
      );
      footer$srcFiles <- c(footer$srcFiles, list(srcFile=srcFile));
      csvList[[key]] <- csv;
      rm(csv);
    }
    rm(tags);
  }
}
print(csvList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import UFL from CSV files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tags <- sprintf("na%s,%s%s", naVersion, user, datestamp);
ufl <- NULL;
tryCatch({
  ufl <- AromaUflFile$byChipType(getChipType(cdf), tags=tags);
}, error = function(ex) {})
if (is.null(ufl)) {
  ufl <- AromaUflFile$allocateFromCdf(cdf, tags=tags, nbrOfEnzymes=nbrOfEnzymes);
}
print(ufl);

for (kk in seq(along=csvList)) {
  csv <- csvList[[kk]];
  print(csv);
  units <- importFrom(ufl, csv, verbose=log);
  str(units);
  ## GenomeWideSNP_6.na26.annot.csv:    int [1:934968] 334945 334944 ...
  ## GenomeWideSNP_6.cn.na26.annot.csv: int [1:945826] 935622 935777 ...
}

footer <- c(readFooter(ufl), footer);
writeFooter(ufl, footer);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Statistics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
print(ufl);

x <- summaryOfUnits(ufl, enzymes=c("NspI", "StyI"));
print(x);
## GenomeWideSNP_6,Full,na26,HB20080722:
##                 snp    cnp affxSnp other   total
## enzyme1-only 246080 451191       0     0  697271
## enzyme2-only 160899      0       0     0  160899
## both         522472 494615       0     0 1017087
## missing        2495     20    3022   621    6158
## total        931946 945826    3022   621 1881415
