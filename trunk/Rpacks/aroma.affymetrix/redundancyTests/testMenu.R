if (interactive())
  savehistory();
library(aroma.affymetrix);

# Use special file cache for testing
## options("R.cache::rootPath"="~/.Rcache,scratch");
## options("R.cache::touchOnLoad"=TRUE);

source("textUI.R");
source("patch.R");

pathname <- textSelectFile("testScripts", history=TRUE);
if (!is.null(pathname)) {
  patchPackage("aroma.affymetrix");
  source(pathname, echo=TRUE);
}
