setConstructorS3("AptProbesetGenotype", function(..., tags=c("APT", "BRLMM")) {
  extend(AptProbesetGenotype(..., tags=tags), "BrlmmModel")
}, private=TRUE)



############################################################################
# HISTORY:
# 2007-03-22
# o Created "wrapper" class for the BRLMM model, because the class 
#   AptProbesetGenotype is more designed to call the APT binaries.
############################################################################
