if (interactive()) savehistory();
library("aroma.cn");

log <- verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

tagsList <- list("MSKCC", "Harvard", "Stanford", "Broad");
dataSet <- "TCGA,GBM,testSet,pairs";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load raw CN data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dsList <- lapply(tagsList, FUN=function(tags) {
  AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType="*");
});
# Keep only common samples (just in case)
names <- Reduce(intersect, lapply(dsList, FUN=getNames));
dsList <- lapply(dsList, FUN=extract, names);
print(dsList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Extract subset of raw CNs across platforms
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Sample of interest
array <- 4;

# Chromosomal region of interest
chromosome <- 3;
region <- c(0,60)*1e6;
region <- NULL

cnList <- lapply(dsList, FUN=function(ds) {
  df <- getFile(ds, array);
  extractRawCopyNumbers(df, chromosome=chromosome, region=region);
})
print(cnList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot data along chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
xlab <- "Physical position (in Mb)";
Mlab <- expression(M==log[2](theta[T]/theta[N]));

# Find the range physical positions.
xlim <- range(sapply(cnList, FUN=xRange));
Mlim <- c(-1,1)*4;

# The name of the tumor/normal pair
name <- getNames(dsList[[1]])[array];

layout(seq(along=cnList));
par(mar=c(4.2,4.2,1.3,2.1));
for (kk in seq(along=cnList)) {
  cn <- cnList[[kk]];
  plot(cn, xlim=xlim, ylim=Mlim);
  stext(side=3, pos=0, cex=0.8, name);
  stext(side=3, pos=1, cex=0.8, sprintf("Chr%02d (n=%d)", chromosome, nbrOfLoci(cn)));
}
