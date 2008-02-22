.patchCdfMergeStrands <- function() {
  # Nothing to do?
  ver <- packageDescription("affxparser")$Version;
  if (compareVersion(ver, "1.11.5") >= 0)
    return(); 

  cdfMergeStrands.patch <- function(groups, ...) {
    nbrOfGroups <- length(groups);
    if (nbrOfGroups != 2 && nbrOfGroups != 4 && nbrOfGroups %% 2 != 0)
      return(groups);
  
    names <- names(groups);
    unames <- unique(names);
  
    res <- list();
  
    # For each allele...
    for (name in unames) {
      idx <- which(name == names);
      group <- .subset2(groups, idx[1]);
      nfields <- length(group);
      if (nfields > 0) {
        ffs <- 1:nfields;
        idx <- idx[-1];
        while(length(idx) > 0) {
          groupNext <- .subset2(groups, idx[1]);
    
          # For each field...
          for (ff in ffs) {
            fields <- .subset2(group, ff);
            fieldsNext <- .subset2(groupNext, ff);
            ndim <- length(dim(fields));
            if (ndim <= 1) {
              fields <- append(fields, fieldsNext);
            } else if (ndim == 2) {
              fields <- cbind(fields, fieldsNext);
            } else {
              # This should never occur for a normal CDF structure.
              fields <- append(fields, fieldsNext);
            }
            group[[ff]] <- fields;
          }
    
          idx <- idx[-1];
        } # while(...)
      }
  
      res[[name]] <- group;
    } # for (name ...)
  
    res;
  } # cdfMergeStrands()
  
  # Patch affxparser::cdfMergeStrands()
  reassignInPackage("cdfMergeStrands", "affxparser", cdfMergeStrands.patch); 
} # .patchCdfMergeStrands()
