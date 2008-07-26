verbose && enter(verbose, "Adding data sets");
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# CRMA: Without probe-sequence normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sets <- addRocDataSet(sets, "ACC,ra,-XY,AVG,+300,A+B,FLN,-XY", name="CRMA*");
sets <- addRocDataSet(sets, "ACC,ra,-XY,RMA,+300,A+B,FLN,-XY", name="CRMA");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# CRMA: With probe-sequence normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sets <- addRocDataSet(sets, "ACC,ra,-XY,BCN,-XY,+300,RMA,A+B,FLN,-XY", name="CRMA(BCN)");
sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,RMA,A+B,FLN,-XY", name="CRMA+");
sets <- addRocDataSet(sets, "ACC,ra,-XY,BPN,-XY,+300,AVG,A+B,FLN,-XY", name="CRMA+*");

### sets <- addRocDataSet(sets, "ACC,ra,-XY,BCN,-XY,AVG,+300,A+B,FLN,-XY", name="CRMAv0");
### sets <- addRocDataSet(sets, "ACC,ra,-XY,BCN,-XY,RMA,+300,A+B,FLN,-XY", name="CRMAv0");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Other methods
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sets <- addRocDataSet(sets, "GTC", name="GTC");
sets <- addRocDataSet(sets, "APT,Full,QN,sketch,PLIER", name="APT");
sets <- addRocDataSet(sets, "dChip,IN,MBEI,A+B", name="dChip");

verbose && exit(verbose);

verbose && enter(verbose, "Keeping data sets of interest");
# Main ROC comparison
sets <- sets[c("CRMA+", "GTC", "APT", "dChip")];

# Smoothed ROC comparison
sets <- sets[c("CRMA+", "GTC")];

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
