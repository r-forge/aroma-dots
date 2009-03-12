if (interactive()) savehistory();
library("R.rsp");
graphics.off();

sampleName <- "TCGA-06-0154";
sampleName <- "TCGA-12-0616";
sampleName <- "TCGA-12-0618";
sampleName <- "TCGA-12-0620";
sampleName <- "TCGA-08-0386";
outPathname <- sprintf("%s,TumorBoost.pdf", sampleName);
if (!isFile(outPathname)) {
  x <- rsptex("main.tex.rsp", force=TRUE);
  x <- "main.pdf";
  file.rename(x, outPathname);
}
