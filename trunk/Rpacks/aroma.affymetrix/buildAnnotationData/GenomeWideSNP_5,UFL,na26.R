if (interactive()) savehistory();
library("aroma.affymetrix");
log <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
naVersion <- "26";
user <- "HB";
datestamp <- "20080723";

chipType <- "GenomeWideSNP_5";
cdfTags <- "Full,r2";
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
      print(csv);
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
  ## GenomeWideSNP_5.na26.annot.csv:    ...
  ## GenomeWideSNP_5.cn.na26.annot.csv: int [1:312384] 513035 516828 ...
}

footer <- c(readFooter(ufl), footer);
writeFooter(ufl, footer);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Statistics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
x <- summaryOfUnits(ufl, enzymes=c("NspI", "StyI"));
print(x);
## GenomeWideSNP_5,Full,na26,HB20080723:
##                 snp    cnp affxSnp other  total
## enzyme1-only 116979 140099       0     0 257078
## enzyme2-only  74135   1208       0     0  75343
## both         248980 171077       0     0 420057
## missing       60474 104885    3022    69 168450
## total        500568 417269    3022    69 920928
