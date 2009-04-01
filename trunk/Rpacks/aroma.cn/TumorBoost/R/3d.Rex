if (interactive()) savehistory();
library("rgl");
library("aroma.affymetrix");
sourceDirectory("R/");

dataSet <- "TCGA_GBM,Broad,2.0";
dataSetTags <- "ACC,ra,-XY,USE.ME=TRUE,BPN,-XY,AVG,FLN,-XY";
chipType <- "GenomeWideSNP_6";

sampleName <- "TCGA-06-0154";
sampleName <- "TCGA-12-0616";
sampleName <- "TCGA-12-0618";
sampleName <- "TCGA-08-0386";
sampleName <- "TCGA-12-0620";
sampleName <- "TCGA-02-0069";
chromosome <- 4;

if (!exists("ces") || !inherits(ces, "SnpChipEffectSet")) {
  ces <- SnpChipEffectSet$byName(dataSet, tags=dataSetTags, 
                                   chipType=chipType, mergeStrands=TRUE);
  setFullNamesTranslator(ces, function(names, ...) {
    pattern <- "^(TCGA-[0-9]{2}-[0-9]{4})-([0-9]{2}.)-(.*)$"
    names <- gsub(pattern, "\\1,\\2,\\3", names);
    names <- gsub(",chipEffects", "", names);
    names;
  })

  print(ces);
}

# Extract all normals 
is10 <- sapply(ces, FUN=function(cf) any(regexpr("^10", getTags(cf)) != -1));
is11 <- sapply(ces, FUN=function(cf) any(regexpr("^11", getTags(cf)) != -1));
isNormal <- (is10 | is11);
cesN <- extract(ces, isNormal);
ceNR <- getAverageFile(cesN);
print(ceNR);

cdf <- getCdf(ces);
unitTypes <- getUnitTypes(cdf);
ugp <- AromaUgpFile$byChipType(getChipType(ces));

# Extract (tumor,normal) pair
cesTN <- extract(ces, indexOf(ces, sampleName));
typeTags <- sapply(cesTN, getTags)[1,];
cesTN <- extract(cesTN, order(typeTags));
stopifnot(length(cesTN) == 2);

# Extract tumor normal pair of interest
pairTag <- sapply(cesTN, FUN=function(df) getTags(df)[1]);
pairTag <- paste(pairTag, collapse="vs");
fullname <- sprintf("%s,%s", sampleName, pairTag);



data <- getUnitsToStudy(cesTN, chromosome, subset=NULL);
xyz <- cbind(x=data$position, y=data$thetaC, z=asTotalFracBSnpData(data$thetaTB)[,2]);
xyz <- xyz[order(xyz[,"x"]),];

open3d();
par3d(windowRect=c(50,50,1000,300));
jj <- 1:15000;
plot3d(xyz[jj,], aspect=c(4,1,1)); 
