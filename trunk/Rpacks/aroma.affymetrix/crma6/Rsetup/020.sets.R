verbose && enter(verbose, "Adding data sets");
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# CRMA: Without probe-sequence normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sets <- addRocDataSet(sets, "ACC,ra,-XY,AVG,+300,A+B,FLN,-XY", name="CRMA*");
sets <- addRocDataSet(sets, "ACC,ra,-XY,RMA,+300,A+B,FLN,-XY", name="CRMA");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# CRMA: With probe-sequence normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### # BCN
### sets <- addRocDataSet(sets, "ACC,ra,-XY,BCN,-XY,+300,RMA,A+B,FLN,-XY", name="CRMA(BCN)");
### sets <- addRocDataSet(sets, "ACC,ra,-XY,BCN,-XY,AVG,+300,A+B,FLN,-XY", name="CRMAv0");
### sets <- addRocDataSet(sets, "ACC,ra,-XY,BCN,-XY,RMA,+300,A+B,FLN,-XY", name="CRMAv0");
### 
### # BPN
### sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,RMA,A+B,FLN,-XY", name="CRMA+");
### sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,df=9,RMA,A+B,FLN,-XY", name="DRMA+");
### sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+200,RMA,A+B,FLN,-XY", name="CRMA+200");
### sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,RMA,A+B,FLN,-XY,UTSN,-XY,v2", name="CRMA#");
### sets <- addRocDataSet(sets, "RBC,BPN,-XY,RMA,A+B,FLN,-XY", name="CRMA-");
###                              
### # BPN (single-array)
### sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,AVG,A+B,FLN,-XY", name="CRMA+*");
### sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,df=9,AVG,FLN,-XY", name="DRMA+*");
### sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,AVG,A+B,FLN,-XY,UTSN,-XY", name="CRMA+*#");



# BPN,z
sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,z,RMA,A+B,FLN,-XY,UTSN,-XY", name="FRMA+");
sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,z,RMA,A+B,FLN,-XY", name="ERMA+");
sets <- addRocDataSet(sets, "BPN,-XY,+300,z,ACC,ra,-XY,RMA,A+B,FLN,-XY,UTSN,-XY", name="CRMA3");
sets <- addRocDataSet(sets, "BPN,-XY,+300,z,ACC,ra,-XY,RMA,A+B,FLN,-XY", name="CRMA4");
sets <- addRocDataSet(sets, "BPN,-XY,z,ACC,ra,-XY,RMA,+300,A+B,FLN,-XY,UTSN,-XY", name="CRMA5");
sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,z,RMA,A+B,FLN,-XY,z", name="GRMA+");

# BPN,z (single-array)
sets <- addRocDataSet(sets, "BPN,-XY,z,ACC,ra,-XY,AVG,FLN,-XY,+300", name="CRMA5*");
sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,z,AVG,FLN,-XY", name="ERMA+*");
sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,z,AVG,A+B,FLN,-XY,z", name="GRMA+*");




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Other methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sets <- addRocDataSet(sets, "GTC", name="GTC");
sets <- addRocDataSet(sets, "APT,Full,QN,sketch,PLIER", name="APT");
sets <- addRocDataSet(sets, "dChip,IN,MBEI,A+B", name="dChip0");
sets <- addRocDataSet(sets, "dChip,MBEI", name="dChip");
sets <- addRocDataSet(sets, "dChip,AVG", name="dChip*");

sets <- sets[!is.na(names(sets))];
verbose && print(verbose, names(sets));
verbose && exit(verbose);

verbose && enter(verbose, "Keeping data sets of interest");
# Main ROC comparison
keep <- c("CRMA+", "GTC", "APT", "dChip");

# Smoothed ROC comparison

#keep <- c("CRMA+*#", "CRMA+*", "CRMA+", "GTC");

keep <- c("CRMA+", "CRMA", "GTC", "dChip");
keep <- c("ERMA+", "CRMA+", "ERMA+*", "CRMA+*", "GTC", "dChip");
keep <- c("ERMA+", "CRMA+", "ERMA+*", "CRMA+*", "GTC", "dChip");
keep <- c("CRMA3", "CRMA4", "FRMA+", "ERMA+", "CRMA+", "ERMA+*", "GTC");
keep <- c("CRMA6", "CRMA5", "ERMA+", "CRMA+", "GTC");
keep <- c("GRMA+", "ERMA+*", "GTC");
keep <- c("GRMA+", "ERMA+*", "GTC", "dChip", "dChip*");
keep <- c("GRMA+", "ERMA+*", "GRMA+*");

hasPrefix <- function(name, prefix, ...) {
  (substring(name, 1, nchar(prefix)) == prefix);
}
keep2 <- sapply(names(sets), FUN=function(name) {
  any(sapply(keep, FUN=function(prefix) hasPrefix(name, prefix)));
});
#sets <- sets[keep2];
sets <- sets[keep];
verbose && print(verbose, names(sets));
verbose && exit(verbose);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Load data and calculate log-ratios, if missing
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
verbose && enter(verbose, "Adding ROC data");
for (kk in seq(along=sets)) {
  key <- names(sets)[kk];
  set <- sets[[key]];
  set <- addRocData(set, verbose=verbose);
  sets[[key]] <- set;
  rm(set);
}
verbose && exit(verbose);


verbose && enter(verbose, "Updating graphics annotations");
sets <- updateGraphics(sets);
verbose && exit(verbose);

nbrOfDataSets <- sum(sapply(sets, FUN=function(set) !is.null(set$pathname)));
