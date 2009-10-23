  # "Extract Name":
  #  Examples:
  #  (1) TCGA-02-0001-01C-01D
  #  (2) TCGA-02-0001-01C-01D-00182-01
  #  Format:
  #  <extract name> := <portion id><a>-<n>*-<nn>, e.g. 'TCGA-02-0001-01C-01D-00182-01'
  #  <portion id> := <sample id>-<nn>, e.g. 'TCGA-02-0001-01C-01'
  #  <sample id> := <patient id>-<nn>, e.g. 'TCGA-02-0001-01C'
  #  <patient id> := <name>-<nn><a>..., e.g. 'TCGA-02-0001'

parseSdrfExtractName <- function(extractName, ...) {
  parts <- unlist(strsplit(extractName, split="-"));

  patientId <- paste(parts[1:3], collapse="-");
  sampleId <- paste(patientId, parts[4], sep="-");
  portionId <- paste(sampleId, gsub("[A-Z]$", "", parts[5]), sep="-");
  centerId <- parts[length(parts)];

  list(extractName=extractName, patientId=patientId, sampleId=sampleId, 
                                    centerId=centerId, portionId=portionId);
} # parseSdrfExtractName()

############################################################################
# HISTORY:
# 2008-03-19
# o Created from EP's code/explainations/my guesses. [HB].
############################################################################

