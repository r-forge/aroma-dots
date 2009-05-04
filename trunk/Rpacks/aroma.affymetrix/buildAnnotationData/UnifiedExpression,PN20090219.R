if (interactive() savehistory();
library("aroma.affymetrix");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# BioMart data file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
chipType <- "UnifiedExpression";
pathname <- findAnnotationDataByChipType(chipType, pattern="biomaRtPositions");
db <- TabularTextFile(pathname);
setColumnNameTranslator(db, function(names, ...) {
  names <- gsub("id", "unitName", names, fixed=TRUE);
  names <- gsub("chromosome_name", "chromosome", names, fixed=TRUE);
  names <- gsub("_", " ", names, fixed=TRUE);
  names <- toCamelCase(names);
  names;
})
print(db);
nbrOfUnits <- nbrOfRows(db);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Read data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
df <- readDataFrame(db, colClassPattern=c("*"="integer", "^unitName$"="character"), na.strings="NA");
str(df);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Allocate UGP file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
tags <- c("geneLevel", "HB20090219");
filename <- sprintf("%s.ugp", paste(c(chipType, tags), collapse=","));
ugp <- AromaUgpFile$allocate(filename=filename, path=getPath(db), nbrOfRows=nbrOfUnits, platform="TCGA-Multisource", chipType=chipType);
footer <- readFooter(ugp);
footer$srcFile <- list(filename=getFilename(db), checksum=getChecksum(db));
footer$author <- list(name="Pierre Neuvial", email="pierre@stat.berkeley.edu");
writeFooter(ugp, footer);
print(ugp);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp[,1] <- df[,"chromosome"];
ugp[,2] <- df[,"position"];


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
gi <- UgpGenomeInformation$byChipType(chipType, tags="geneLevel,HB20090219");
print(gi);
print(getChromosomeStats(gi));
