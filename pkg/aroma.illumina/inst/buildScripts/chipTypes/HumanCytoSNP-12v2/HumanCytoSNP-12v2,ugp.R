if (interactive()) savehistory();
library("aroma.core");
verbose <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
platform <- "Illumina";
chipType <- "HumanCytoSNP-12v2";
tags <- "HB20110924";

createdBy <- "Henrik Bengtsson (hb@biostat.ucsf.edu)";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Unit names mapping
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
unf <- TextUnitNamesFile$byChipType(chipType, tags=tags);
print(unf);

unitNames <- getUnitNames(unf);
printf("Unit names:\n");
str(unitNames);


filename <- "marker_list_humancytoSNP12.txt";
db <- TabularTextFile(filename, path=getPath(unf));
print(db);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Allocate UNF file (gives an error if already exists)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp <- AromaUgpFile$allocateFromUnitNamesFile(unf, tags=tags);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import unit names
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

patterns <- c("*"="NULL", Name="character", Chr="character", Position="integer");
data <- readDataFrame(db, colClassPatterns=patterns);

nbrOfUnits <- nrow(data);
printf("Number of units: %d\n", nbrOfUnits);

# Translate unit names to unit indices via UNF
units <- indexOf(unf, names=data[,"Name"]);
printf("Units:\n");
str(units);

# Sanity check
units <- Arguments$getIndices(units, max=nbrOfUnits(ugp));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Translate chromosome ids to chromosome indices
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
map <- c(1:24);
names(map) <- map;
map <- c(map, X=23, Y=24, M=25, XY=26);
data[,"Chr"] <- map[data[,"Chr"]];
printf("Number of units per chromosome:\n");
print(table(data[,"Chr"]));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Write annotation data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp[,1] <- data[,"Chr"];
ugp[,2] <- data[,"Position"];

ftr <- readFooter(ugp);
ftr$createdBy <- createdBy
ftr$srcFile <- list(filename=getFilename(db), filesize=getFileSize(db), checksum=getChecksum(db));
writeFooter(ugp, ftr);

print(ugp);
