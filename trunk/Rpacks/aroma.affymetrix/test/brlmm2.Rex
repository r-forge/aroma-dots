library(oligo);
library(aroma.affymetrix);
source("init.R");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
.Machine$float.eps <- sqrt(.Machine$double.eps);
tolerance <- .Machine$float.eps^0.5;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Select data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
path <- "raw/Test100K/Mapping50K_Hind240";  # Works! 2006-12-13
path <- "raw/Test500K/Mapping250K_Nsp";

path <- "raw/SinclairA_etal_2006/Mapping250K_Sty";

ds <- AffymetrixCelSet$fromFiles(path);
print(ds);

apt <- AptProbesetGenotype(ds);
callSet <- process(apt, verbose=verbose);
