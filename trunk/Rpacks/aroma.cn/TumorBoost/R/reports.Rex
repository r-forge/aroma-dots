if (interactive()) savehistory();
library("R.rsp");

sampleName <- "TCGA-06-0154";
sampleName <- "TCGA-12-0616";
sampleName <- "TCGA-12-0618";
outPathname <- sprintf("%s,TumorBoost.pdf", sampleName);
if (!isFile(outPathname)) {
  x <- rsptex("main.tex.rsp", force=TRUE);
  x <- "main.pdf";
  file.rename(x, outPathname);
}
