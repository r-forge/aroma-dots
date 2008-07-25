library("R.utils");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

sourceDirectory("Rincl/", modifiedOnly=TRUE);

# Setup of data etc
sourceDirectory("Rsetup/", modifiedOnly=TRUE);

mkdirs(figPath);
