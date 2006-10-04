setMethodS3("fitCrlmm", "SnpChipEffectSet", function(this, minLLRforCalls=c(50,40,50), ..., verbose=FALSE) {
  recalibrate <- TRUE;
  returnCorrectedM <- TRUE;
  returnParams <- TRUE;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose2 <- as.logical(verbose);

  verbose && enter(verbose, "Fitting CRLMM");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract chip effects as a SnpQSet object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting chip effects as a SnpQSet object");
  qs <- extractSnpQSet(ces, verbose=less(verbose));
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Dimension of data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSnps <- dim(qs)[1];
  nbrOfSamples <- dim(qs)[2];
  verbose && printf(verbose, "Number of SNPs: %d", nbrOfSnps);
  verbose && printf(verbose, "Number of samples: %d", nbrOfSamples);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get M corrections
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculate M corrections using fitAffySnpMixture()");
  correction <- fitAffySnpMixture(qs, verbose=verbose2);
  verbose && str(verbose, correction);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get gender
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving genders");
  if(is.null(qs$gender)){
    verbose && enter(verbose, "Predicting genders from data");
    qs$gender <- snpGenderCall(qs);
    verbose && exit(verbose);
  }
  maleIndex <- (qs$gender == "male");
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cleanCdfName <- annotation(qs);
  pkgName <- sprintf("pd%s", cleanCdfName);
  if (identical(packageDescription(pkgName), NA))
    throw("Package not installed: ", pkgName);

  verbose && enter(verbose, "Loading a priori CRLMM parameters");
  verbose && cat(verbose, "Source package: ", pkgName);
  filename <- sprintf("%sCrlmmInfo.rda", cleanCdfName);
  pathname <- system.file(file.path("data", filename), package=pkgName);
  verbose && cat(verbose, "Pathname: ", pathname);
  verbose && printf(verbose, "File size: %.2fMb", 
                                           file.info(pathname)$size/1024^2);
  vars <- load(pathname);
  verbose && cat(verbose, "Loaded variable(s): ", paste(vars, collapse=", "));
  verbose && exit(verbose);

  # Get CRLMM parameters enviroment
  var <- sprintf("%sCrlmm", cleanCdfName);
  if (!var %in% vars) {
    throw("Internal error. Variable '", var, "' was not among the loaded: ",
          paste(vars, collapse=", "));
  }
  env <- get(var, mode="environment");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the indices of SNPs that are not called by the HapMap project 
  # or that were poorly called.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  index <- which(!env$hapmapCallIndex | env$badCallIndex);
  verbose && cat(verbose, "Among the ", nbrOfSnps, " SNPs, there are ", length(index), " SNPs that were not called by the HapMap project or were poorly called by it.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get initial calls for non-HapMap calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting initial calls for non-HapMap calls");
  myCalls <- matrix(NA, nrow=nbrOfSnps, ncol=nbrOfSamples);
  myCalls[index,] <- getInitialAffySnpCalls(correction, index, 
                                                         verbose=verbose2);
  verbose && exit(verbose);
  gc();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get SNP genotype region parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting SNP genotype region parameters");
  rparams <- getAffySnpGenotypeRegionParams(qs, myCalls,
                             correction$fs, subset=index, verbose=verbose2);
  # Not needed anymore
  rm(myCalls); gc();
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Shrink genotype regions toward the priors
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Shrinking genotype regions toward the priors");
  rparams <- updateAffySnpParams(rparams, env$priors, verbose=verbose2);
  verbose && exit(verbose);
  gc();
  
  verbose && enter(verbose, "Replacing region parameters");
  params <- replaceAffySnpParams(env$params, rparams, index)
  verbose && exit(verbose);
  gc();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate call distances
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating call distances");
  myDist <- getAffySnpDistance(qs, params, correction$fs);
  verbose && exit(verbose);

  XIndex <- getChrXIndex(qs);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the calls from the new distances
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Update the calls from the new distances");
  myCalls <- getAffySnpCalls(myDist, XIndex, maleIndex, verbose=verbose2);
  gc();
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate log-likehood ratios
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating log-likehood ratios");
  LLR <- getAffySnpConfidence(myDist, myCalls, XIndex, maleIndex, 
                                                         verbose=verbose2);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Recalibrate?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (recalibrate) {
    verbose && enter(verbose, "Recalibrating");

    rm(myDist, rparams); gc(); # Not needed anymore

    verbose && enter(verbose, "Resetting poor genotype calls");
    # For each SNP
    for (ii in 1:nrow(myCalls)) {
      # For each genotype class
      for (kk in 1:3) {
        bad <- which(LLR[ii,] < minLLRforCalls[kk] & myCalls[ii,] == kk);
        myCalls[ii, bad] <- NA;
      }
      rm(bad);
    }
    rm(LLR); gc(); # Not needed anymore 
    verbose && exit(verbose);

    verbose && enter(verbose, "Getting SNP genotype region parameters");
    rparams <- getAffySnpGenotypeRegionParams(qs, myCalls, correction$fs,
                                                         verbose=verbose2);
    rm(myCalls); gc(); # Not needed anymore
    verbose && exit(verbose);

    verbose && enter(verbose, "Update the calls from the new distances");
    rparams <- updateAffySnpParams(rparams, env$priors);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating call distances");
    myDist <- getAffySnpDistance(qs, rparams, correction$fs,
                                                          verbose=verbose2);
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Update the calls from the new distances");
    myCalls <- getAffySnpCalls(myDist, XIndex, maleIndex, verbose=verbose2);
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Calculating log-likehood ratios");
    LLR <- getAffySnpConfidence(myDist, myCalls, XIndex, maleIndex,
                                                         verbose=verbose2);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # if (recalibrate)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create return structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ret <- list(calls=myCalls, llr=LLR);
  rm(LLR); gc(); # Not needed anymore


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Corrected log-ratios
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (returnCorrectedM) {
    SIGN <- c(+1,0,-1);
    cM <- getM(qs);
    for(ii in 1:nrow(cM)){
      # For the two homozygote groups
      for(kk in c(1,3)) {
        index <- which(myCalls[ii,] == kk);
        cM[ii,index,] <- cM[ii,index,] + 
                             SIGN[kk]*(rparams$f0-correction$fs[ii,kk,]);
      }
    }
    ret$M <- cM;
    rm(cM, index, SIGN); gc(); # Not needed anymore
  } # if (returnCorrectedM)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Feature data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (returnParams) {
    rparams <- as.data.frame(rparams);
    fD <- new("AnnotatedDataFrame",
      data=rparams,
      varMetadata=data.frame(labelDescription=colnames(rparams),
                                                row.names=colnames(rparams))
    )
  } else {
    fD <- new("AnnotatedDataFrame");
  }
  
  gender <- rep("female", length(maleIndex));
  gender[maleIndex] <- "male";
  rm(maleIndex); gc(); # Not needed anymore

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Phenotype data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(qs$gender)) {
    phenoData <- new("AnnotatedDataFrame",
                        data=cbind(pData(qs),
                          data.frame(crlmmSNR=as.numeric(correction$snr),
                                     gender=gender,
                                     row.names=sampleNames(qs))),
                        varMetadata= rbind(varMetadata(qs),
                          data.frame(labelDescription=c("crlmmSNR", "gender"),
                                     row.names=c("crlmmSNR", "gender"))))
  } else {
    phenoData <- new("AnnotatedDataFrame",
                        data=cbind(pData(qs),
                          data.frame(crlmmSNR=as.numeric(correction$snr),
                                     row.names=sampleNames(qs))),
                        varMetadata= rbind(varMetadata(qs),
                          data.frame(labelDescription=c("crlmmSNR"),
                                     row.names=c("crlmmSNR"))))
  }
  rm(correction); gc(); # Not needed anymore


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return genotype estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- new("SnpCallSetPlus",
              featureData=fD,
              phenoData=phenoData,
              experimentData=experimentData(qs),
              annotation=annotation(qs),
              calls=ret$calls,
              callsConfidence=ret$llr,
              logRatioAntisense=ret$M[,,"antisense"],
              logRatioSense=ret$M[,,"sense"]);
  verbose && exit(verbose);

  res;
}) # fitCrlmm()


############################################################################
# HISTORY:
# 2006-10-03
# o Created based on oligo::crlmm() v0.99.20.
############################################################################
