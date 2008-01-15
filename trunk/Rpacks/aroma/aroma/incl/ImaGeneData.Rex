######################################################################
# Example that reads two ImaGene data files, which are single-channel
# and represent one slide, and combines them into a RawData object.
######################################################################
igG <- ImaGeneData$read("imagene234.cy3", path=system.file("data-ex", package="aroma"))
igR <- ImaGeneData$read("imagene234.cy5", path=system.file("data-ex", package="aroma"))

# Get the raw data; Rfg, Gfg, Rbg and Gbg
raw <- getRawData(igR, igG)

# The the background corrected data
ma <- getSignal(raw)

# Plot M vs A with a lowess line through the data points
subplots(4)
plot(ma)
plotSpatial(ma)
boxplot(ma, groupBy="printtip")


######################################################################
# Example that reads multiple ImaGene data files, which all are
# single-channel and represent several slide, and then combines them
# into a RawData object of several slides.
######################################################################

# Function that reads several (single-channel) ImaGene data files
# and combines them into a multi-slide two-channel RawData object
myImaGeneReadAll <- function(slides, extensions=c("cy3","cy5"), ...) {
  # An empty RawData object to store all slides
  raw <- RawData()
  for (slide in slides) {
    igCy3 <- ImaGeneData$read(paste(slide, extensions[1], sep="."), ...)
    igCy5 <- ImaGeneData$read(paste(slide, extensions[2], sep="."), ...)
    tmp <- getRawData(igCy3, igCy5)
    append(raw, tmp)
  }

  # Return the RawData object containing all slides
  raw;
} # myImaGeneReadAll()


# "Faking" three slides...
slides <- rep("imagene234", 3)

# ...and reads them.
raw <- myImaGeneReadAll(slides, path=system.file("data-ex", package="aroma"))
