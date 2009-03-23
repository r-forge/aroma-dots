if (interactive() savehistory();
library("aroma.affymetrix");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Annotation data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
cdf <- AffymetrixCdfFile$byChipType("HG-U133_Plus_2");
print(cdf);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Allocate UGP file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp <- AromaUgpFile$allocateFromCdf(cdf, tags=("center,MO20090219"));
print(ugp);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Data to import
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
path <- file.path(getPath(cdf), "OrtizM");
filename <- "Human Genome U133 Plus 2.0 Array genome info.xls";
db <- TabularTextFile(filename, path=path);
print(db);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Read and tranform data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
df <- readDataFrame(db, colClassPattern=c("*"="character", "^Pos"="double"));
str(df);

print(table(df$chromosome));
map <- c(X=23, Y=24, M=25);
for (kk in seq(along=map)) {
  df[,"chromosome"] <- gsub(names(map)[kk], map[kk], df[,"chromosome"]);
}
df[,"chromosome"] <- as.integer(df[,"chromosome"]);
print(table(df[,"chromosome"]));

# Map rows in table to unit indices
units <- indexOf(cdf, names=df[,"Probe Set"]);
stopifnot(all(is.finite(units)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Store UGP data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ugp[units,1] <- df[,"chromosome"];
ugp[units,2] <- df[,"Position"];


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update file footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
footer <- readFooter(ugp);
footer$srcFile <- list(filename=getFilename(db), checksum=getChecksum(db));
footer$createdBy <- list(name="Henrik Bengtsson", email="hb@stat.berkeley.edu", comment="Recieved file Human Genome U133 Plus 2.0 Array genome info.xls from Maria Ortiz, CEIT, Spain");
writeFooter(ugp, footer);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
gi <- getGenomeInformation(cdf);
print(gi);
print(getChromosomeStats(gi));
