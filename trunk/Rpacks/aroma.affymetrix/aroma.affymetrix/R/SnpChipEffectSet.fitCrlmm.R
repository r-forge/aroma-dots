###########################################################################/**
# @set "class=SnpChipEffectSet"
# @RdocMethod fitCrlmm
#
# @title "Fits chip effects to the CRLMM model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{minLLRforCalls}{A @numeric @vector of length three.}
#   \item{recalibrate}{If @TRUE, a second round of the fitting is done.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "oligo::SnpCallSetPlus-class" object.
# }
#
# \details{
#   Simple benchmarking on shadowfax: 
#   For 90 CEPH Xba chips it takes ~90 minutes, i.e. 60 s/array.
# }
#
# @author
#
# \seealso{
#   @see "oligo::crlmm".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fitCrlmm", "SnpChipEffectSet", function(this, minLLRforCalls=c(AA=50, AB=40, BB=50), recalibrate=TRUE, transform=c("log", "asinh"), ..., verbose=FALSE) {
  returnCorrectedM <- TRUE;
  returnParams <- TRUE;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'minLLRforCalls':
  minLLRforCalls <- Arguments$getDoubles(minLLRforCalls, length=3);

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
  qs <- extractSnpQSet(this, transform=transform, verbose=less(verbose));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Dimension of data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSnps <- dim(qs)[1];
  nbrOfSamples <- dim(qs)[2];
  verbose && printf(verbose, "Number of SNPs: %d\n", nbrOfSnps);
  verbose && printf(verbose, "Number of samples: %d\n", nbrOfSamples);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get gender
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # getChrXIndex() is robust against SNP reordering /HB 2006-10-03
  xIndex <- getChrXIndex(qs);
  
  verbose && enter(verbose, "Retrieving genders");
  if(is.null(qs$gender)){
    verbose && enter(verbose, "Predicting genders from data");
    # snpGenderCall() is robust against SNP reordering /HB 2006-10-03
    qs$gender <- snpGenderCall(qs, xIndex=xIndex);
    verbose && exit(verbose);
  }
  maleIndex <- (qs$gender == "male");
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get M corrections
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving M corrections");
  # It *looks like* this is robust against SNP reordering /HB 2006-10-03
  key <- list(method="fitAffySnpMixture", path=getPath(this), 
                                         sampleNames=getSampleNames(this));
  correction <- loadCache(key=key);
  if (is.null(correction)) {
    verbose && enter(verbose, "Calculating M corrections using fitAffySnpMixture()");
    correction <- fitAffySnpMixture(qs, xIndex=xIndex, verbose=verbose2);
    verbose && exit(verbose);
    comment <- paste(unlist(key, use.names=FALSE), collapse=";");
    saveCache(key=key, correction, comment=comment);
  }
  verbose && str(verbose, correction);
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
  verbose && printf(verbose, "File size: %.2fMB\n", 
                                           file.info(pathname)$size/1024^2);
  vars <- load(pathname);
  verbose && cat(verbose, "Loaded variable(s): ", paste(vars, collapse=", "));

  # Get CRLMM parameters enviroment
  var <- sprintf("%sCrlmm", cleanCdfName);
  verbose && enter(verbose, "Getting variable '", var, "'");
  if (!var %in% vars) {
    throw("Internal error. Variable '", var, "' was not among the loaded: ",
          paste(vars, collapse=", "));
  }
  env <- get(var, mode="environment");
  verbose && exit(verbose);
  rm(list=vars); # Not needed anymore
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the indices of SNPs that are not called by the HapMap project 
  # or that were poorly called.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying poor or missing HapMap calls");
  snpNames <- featureNames(qs);
  verbose && cat(verbose, "Feature names:");
  verbose && str(verbose, snpNames);
  value <- !env$hapmapCallIndex;
  poorCalls <- value[match(snpNames, names(value))];
  value <- env$badCallIndex;
  poorCalls <- poorCalls | value[match(snpNames, names(value))];
  rm(value);
  poorCalls <- which(poorCalls);
  verbose && cat(verbose, "Among the ", nbrOfSnps, " SNPs, there are ", length(poorCalls), " SNPs that were not called by the HapMap project or were poorly called by it.");
  verbose && exit(verbose);

  # Get rid of large enviroment again
  priors <- env$priors;
  params <- env$params;
  rm(env); gc(); # Not needed anymore

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get initial calls for non-HapMap calls
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting initial calls for non-HapMap calls");
  myCalls <- matrix(NA, nrow=nbrOfSnps, ncol=nbrOfSamples);
  myCalls[poorCalls,] <- getInitialAffySnpCalls(correction, 
                                       subset=poorCalls, verbose=verbose2);
  verbose && exit(verbose);
  gc();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get SNP genotype region parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting SNP genotype region parameters");
  rparams <- getAffySnpGenotypeRegionParams(qs, myCalls,
                        correction$fs, subset=poorCalls, verbose=verbose2);
  # Not needed anymore
  rm(myCalls); gc();
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Shrink genotype regions toward the priors
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Shrinking genotype regions toward the priors");
  rparams <- updateAffySnpParams(rparams, priors, verbose=verbose2);
  params <- replaceAffySnpParams(params, rparams, subset=poorCalls);
  verbose && exit(verbose);
  gc();
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate call distances
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating call distances");
  # SNPs in qs, params, and correction$fs must be ordered the same.
  # myDist is a JxIx3x2 matrix where J=#SNPs, I=#arrays, with
  # 3=#genotypes and # 2=#strands.
  myDist <- getAffySnpDistance(qs, params, correction$fs);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the calls from the new distances
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Updating the calls from the new distances");
  # This takes time: This is worth optimizing.
  res <- list(myDist=myDist, xIndex=xIndex, maleIndex=maleIndex);
  saveCache(key=list("yo"), res);
  myCalls <- getAffySnpCalls(myDist, xIndex, maleIndex, verbose=verbose2);
  gc();
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate log-likehood ratios
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating log-likehood ratios");
  llr <- getAffySnpConfidence(myDist, myCalls, xIndex, maleIndex, 
                                                         verbose=verbose2);
  verbose && str(verbose, llr);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Recalibrate?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (recalibrate) {
    verbose && enter(verbose, "Recalibrating");

    rm(myDist, rparams); gc(); # Not needed anymore

    verbose && enter(verbose, "Resetting poor genotype calls");
    # For each genotype class
    for (kk in 1:3) {
      bad <- which(myCalls == kk & llr < minLLRforCalls[kk]);
      myCalls[bad] <- NA;
      rm(bad);
    }
    verbose && cat(verbose, "Number of poor calls: ", sum(is.na(myCalls)),
                                              " out of ", length(myCalls));
    rm(llr); gc(); # Not needed anymore
    verbose && exit(verbose);

    verbose && enter(verbose, "Getting SNP genotype region parameters");
    rparams <- getAffySnpGenotypeRegionParams(qs, myCalls, correction$fs,
                                                         verbose=verbose2);
    rm(myCalls); gc(); # Not needed anymore
    verbose && exit(verbose);

    verbose && enter(verbose, "Updating the calls from the new distances");
    rparams <- updateAffySnpParams(rparams, priors);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calculating call distances");
    myDist <- getAffySnpDistance(qs, rparams, correction$fs,
                                                          verbose=verbose2);
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Updating the calls from the new distances");
    myCalls <- getAffySnpCalls(myDist, xIndex, maleIndex, verbose=verbose2);
    verbose && exit(verbose);
    
    verbose && enter(verbose, "Calculating log-likehood ratios");
    llr <- getAffySnpConfidence(myDist, myCalls, xIndex, maleIndex,
                                                         verbose=verbose2);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # if (recalibrate)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create return structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating return structure");
  ret <- list(calls=myCalls, llr=llr);
  rm(myCalls, llr); gc(); # Not needed anymore
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Corrected log-ratios
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (returnCorrectedM) {
    verbose && enter(verbose, "Correcting log-ratios");
    SIGN <- c(+1,0,-1);
    cM <- getM(qs);
    for(ii in 1:nrow(cM)){
      # For the two homozygote groups
      for(kk in c(1,3)) {
        index <- which(ret$calls[ii,] == kk);
        cM[ii,index,] <- cM[ii,index,] + 
                             SIGN[kk]*(rparams$f0-correction$fs[ii,kk,]);
      }
    }
    ret$M <- cM;
    rm(cM, index, SIGN); gc(); # Not needed anymore
    verbose && exit(verbose);
  } # if (returnCorrectedM)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Feature data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (returnParams) {
    rparams <- as.data.frame(rparams);
    featureData <- new("AnnotatedDataFrame",
      data=rparams,
      varMetadata=data.frame(labelDescription=colnames(rparams),
                                                row.names=colnames(rparams))
    )
  } else {
    featureData <- new("AnnotatedDataFrame");
  }
  rm(rparams); gc(); # Not needed anymore


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Phenotype data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating phenotype result");
  gender <- rep("female", nbrOfSamples);
  gender[maleIndex] <- "male";
  rm(maleIndex); gc(); # Not needed anymore
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
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return genotype estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating SnpCallSetPlus object to return");
  res <- new("SnpCallSetPlus",
           featureData=featureData,
           phenoData=phenoData,
           experimentData=experimentData(qs),
           annotation=annotation(qs),
           calls=ret$calls,
           callsConfidence=ret$llr,
           logRatioAntisense=ret$M[,,"antisense"],
           logRatioSense=ret$M[,,"sense"]
         );
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}) # fitCrlmm()


############################################################################
# HISTORY:
# 2006-10-03
# o Created based on oligo::crlmm() v0.99.20.
############################################################################
