library(aroma.affymetrix);
source("init.R");

pngDev <- System$findGraphicsDevice();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the family of datasets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dataSetName <- "CarmichaelC_etal_2006-500k";
#dataSetName <- "Affymetrix_2006-HapMap270";
dataSetName <- "SlaterH_etal_2006";
chipTypes <- c("Mapping250K_Nsp", "Mapping250K_Sty");

#dataSetName <- "SlaterH_etal_2004";
#chipTypes <- c("Mapping50K_Hind240", "Mapping50K_Xba240");

if (!exists("gas", mode="list"))
  gas <- list();

for (chipType in chipTypes) {
  ga <- gas[[chipType]];
  if (is.null(ga)) {
    path <- filePath("raw", dataSetName, chipType);
    print(path);
    ds <- AffymetrixCelSet$fromFiles(path);
    print(ds);
  
#    ga <- RmaCnAnalyzer(ds); 
#    preprocess(ga, subversionTag="QN", verbose=verbose); 

    ga <- MbeiCnAnalyzer(ds); 
    preprocess(ga, subversionTag="ACT.QN", verbose=verbose); 

    setup(ga, verbose=verbose);

    gas[[chipType]] <- ga;
  }

  plotGlad(ga, verbose=verbose);
} # for (chipType ...)

