if (interactive()) savehistory();
library("R.rsp");
graphics.off();
force <- TRUE;

sampleName <- "TCGA-06-0154";
sampleName <- "TCGA-12-0616";
sampleName <- "TCGA-12-0618";
sampleName <- "TCGA-08-0386";
sampleName <- "TCGA-02-0069";
sampleName <- "TCGA-12-0620";
outPathname <- sprintf("%s,TumorBoost.pdf", sampleName);
if (force || !isFile(outPathname)) {
  x <- rsptex("main.tex.rsp", force=TRUE);
  x <- "main.pdf";
  file.rename(x, outPathname);
#  x <- rsptex("egfr.tex.rsp", force=TRUE);
}

