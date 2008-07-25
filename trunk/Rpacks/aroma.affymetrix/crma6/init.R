library("R.utils");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

sourceDirectory("Rincl/", modifiedOnly=TRUE);

# Setup of data etc
setsExists <- (exists("sets", mode="list") && length(sets) > 0);
sourceDirectory("Rsetup/", modifiedOnly=TRUE && setsExists);

mkdirs(figPath);
