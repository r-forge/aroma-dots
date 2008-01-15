library(aroma)

source("MAData.R")
source("Layout.R")
source("Replicates.R")


if (!exists("gpr")) {
  ma <- MAData();
  for (file in c("slide1.gpr", "slide2.gpr", "slide3.gpr", "slide4.gpr")) {
    cat("Loading ", file, "...", sep="");
    gpr <- GenePixData$read(file, path=system.file("data-ex", package="aroma"));
    raw <- getRawData(gpr);
    ma0 <- getSignal(raw);
    ma$append(ma0);
    rm(gpr);
    rm(raw);
    rm(ma0);
    cat("ok\n");
  }
  layout <- getLayout(ma);
}

if (!exists("dups")) {
  dups <- Replicates$fromType("top-bottom", layout);

  # Set the colors of the genes to loop over the following colors:
  colors <- c("red", "green", "blue", "yellow", "purple", "black");
  dups$setParameter("col", colors);

  # Set the symbols of the genes to loop over the following symbols:
  symbols <- c(17, 18, 19, 21, 22, 23, 24);
  dups$setParameter("pch", symbols);
}

ma$getLayout()$setReplicates(dups)
genes <- 1:10;

plot(ma, pch=".", col="black");
highlight(dups, genes);

