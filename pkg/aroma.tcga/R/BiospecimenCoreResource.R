# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# TCGA data set and file name regular expressions
# [1] http://stackoverflow.com/questions/1313934/how-are-nested-capturing-groups-numbered-in-regular-expressions
#
#  patterns <- BiospecimenCoreResource$getBarcodePatterns();
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Biospecimen Core Resource (BCR)
setConstructorS3("BiospecimenCoreResource", function(...) {
  extend(Object(), "BiospecimenCoreResource");
})

setMethodS3("getBarcodePatterns", "BiospecimenCoreResource", function(static, ...) {
  patterns <- list();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Biospecimen Core Resource (BCR)
  # BCR Barcode Patterns
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # ProjectName: Currently only 'TCGA'
  patterns$projectName <- "TCGA";
  # BCR sample collection center code: 00-99
  patterns$siteId <- "[0-9]{2}";
  # Patient ID: 0000-9999
  patterns$patientId <- "[0-9]{4}";
  # Sample type: 01-09 for tumor, 10-19 for normal, 20-29 for contols
  patterns$sampleType <- "[0-9]{2}";
  # Vial count per patient-sample: A, B, C, ...
  patterns$vial <- "[A-Z]";
  # Sample ID: 
  patterns$sampleId <- with(patterns, sprintf("(%s)(%s)", sampleType, vial));
  # Portion ID: Solid tumor are divided into a sequence 100-120 mg sections called portions.
  patterns$portionNbr <- "[0-9]{2}";
  # Analyte code: D=DNA, R=RNA, T=Total RNA, W=WGA Qiagen, G=WGA GenomePlex
  patterns$analyteCode <- "[DRTWG]";
  patterns$portionId <- with(patterns, sprintf("(%s)(%s)", portionNbr, analyteCode));

  # ProjectName-SiteID-PatientID:
  patterns$patient <- with(patterns, sprintf("(%s)-(%s)-(%s)", 
                                             projectName, siteId, patientId));
  
  # TCGA BCR analyte barcode has form:
  #   ProjectName-SiteID-PatientID-SampleID-PortionID
  patterns$analyteBarcode <- with(patterns, sprintf("(%s)-(%s)-(%s)", 
                                            patient, sampleId, portionId));
  
  
  # PlateID ranges from 0001, 0002, ..., 9999 (up to 9999 96 well plates)
  patterns$plateId <- "[0-9]{4}";
  # CenterID defines the TCGA CGCCs or GSCs and ranges from 01-99
  patterns$centerId <- "[0-9]{2}";
  
  # BCR plate barcode: The BCR provides a unique plate barcode for each 96-well 
  # plate delivered to a center. The plate barcode uses the following convention:
  # PlateID-CenterID
  patterns$plateBarcode <- with(patterns, sprintf("(%s)-(%s)", plateId, centerId));
  
  # The BCR plate barcode persists with each cente's data as part of the aliquot 
  # barcode by concatenating the plate barcode with the BCR analyte barcode and 
  # using a dash (-) to separate the two IDs. The compound ID forms the 
  # BCR Aliquot Barcode.
  patterns$aliqoutBarcode <- with(patterns, sprintf("(%s)-(%s)", analyteBarcode, plateBarcode));

  patterns;
}, static=TRUE)


setMethodS3("getBarcodePattern", "BiospecimenCoreResource", function(static, name, ...) {
  patterns <- getBarcodePatterns(static, ...);
  patterns[[name]];  
}, static=TRUE)



############################################################################
# HISTORY:
# 2009-10-02
# o Created.
############################################################################
