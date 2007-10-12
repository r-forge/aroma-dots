if (interactive())
  savehistory();
library(aroma.affymetrix);
source("textUI.R");

#source("init.R");


pathname <- textSelectFile("testScripts");

if (!is.null(pathname))
  source(pathname, echo=TRUE);	
