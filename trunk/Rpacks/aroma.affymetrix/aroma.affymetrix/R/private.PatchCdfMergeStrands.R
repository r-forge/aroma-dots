.patchCdfMergeStrands <- function() {
  # Nothing to do?
  ver <- packageDescription("affxparser")$Version;
  if (compareVersion(ver, "1.11.5") >= 0)
    return(); 

  # Patch affxparser::cdfMergeStrands()
  reassignInPackage("cdfMergeStrands", "affxparser", cdfMergeStrands); 
} # .patchCdfMergeStrands()
