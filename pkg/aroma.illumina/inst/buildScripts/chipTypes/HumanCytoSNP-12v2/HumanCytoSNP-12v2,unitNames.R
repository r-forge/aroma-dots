if (interactive()) savehistory();
library("aroma.core");
verbose <- Verbose(threshold=-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
platform <- "Illumina";
chipType <- "HumanCytoSNP-12v2";

createdBy <- "Henrik Bengtsson <hb@biostat.ucsf.edu>";
tags <- "HB20110924";

fullname <- paste(c(chipType, tags), collapse=",");
filename <- sprintf("%s,unitNames.txt", fullname);
path <- file.path("annotationData", "chipTypes", chipType);
pathname <- Arguments$getWritablePathname(filename, path=path);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Import unit names
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
filename <- "marker_list_humancytoSNP12.txt";
db <- TabularTextFile(filename, path=path);
print(db);

patterns <- c("*"="NULL", Name="character", Chr="character", Position="integer");
patterns <- patterns[1:2];
data <- readDataFrame(db, colClassPatterns=patterns);

nbrOfUnits <- nrow(data);
printf("Number of units: %d\n", nbrOfUnits);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup data frame to be writting
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
hdr <- list(
  chipType = chipType,
  platform = platform,
  nbrOfUnits = nrow(data),
  srcFile = getFilename(db),
  srcFileSize = getFileSize(db),
  srcFileChecksum = getChecksum(db),
  createdBy = createdBy 
);
data <- data.frame(unitName=data[,"Name"]);
writeDataFrame(data, file=pathname, header=hdr);