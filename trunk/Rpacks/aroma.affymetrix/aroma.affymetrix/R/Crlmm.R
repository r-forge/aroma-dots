###########################################################################/**
# @RdocClass Crlmm
#
# @title "The Crlmm class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{ces}{A @see "SnpChipEffectSet" to which this model should 
#     be fitted.}
#   \item{minLLRforCalls}{A @numeric @vector of length three.}
#   \item{transform}{The transform used for the chip effects.}
#   \item{tags}{A @character @vector of tags to be appended to the tags of
#      the input data set.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \seealso{
# }
#*/###########################################################################
setConstructorS3("Crlmm", function(ces=NULL, minLLRforCalls=c(AA=50, AB=40, BB=50), transform=c("log", "asinh"), tags="", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ces':
  if (!is.null(ces)) {
    if (!inherits(ces, "SnpChipEffectSet"))
      throw("Argument 'ces' is not an SnpChipEffectSet object: ", class(ces));
  }

  # Arguments 'minLLRforCalls':
  minLLRforCalls <- Arguments$getDoubles(minLLRforCalls, length=3);

  # Argument 'transform':
  transform <- match.arg(transform);

  # Argument 'tags':
  if (!is.null(ces)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
  }


  extend(Object(), "Crlmm",
    ces = ces,
    minLLRforCalls = minLLRforCalls,
    transform = transform,
    .tags = tags
  )
})



setMethodS3("as.character", "Crlmm", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("SNP chip-effect set: %s", getName(this)));
  ces <- getChipEffects(this);
  tags <- paste(getTags(ces), collapse=",");
  s <- c(s, sprintf("Input tags: %s", tags));
  s <- c(s, sprintf("Transform: %s", this$transform));
  s <- c(s, sprintf("Output tags: %s", paste(getTags(this), collapse=",")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(getCdf(this))));
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
})


###########################################################################/**
# @RdocMethod getRootPath
#
# @title "Gets the root path of this model"
#
# \description{
#  @get "title", which is the class name of the model preceeded by the
#  string \code{"model"}.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getRootPath", "Crlmm", function(this, ...) {
  # Default root path
  paste("model", class(this)[1], sep="");
})



###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the output summarized data set"
#
# \description{
#  @get "title", which is the same as the input data set.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "Crlmm", function(this, ...) {
  ds <- getChipEffects(this);
  getName(ds);
})


###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the output data set"
#
# \description{
#  @get "title", which equals the tags of the input data set plus the tags
#  of this model.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "Crlmm", function(this, ...) {
  ces <- getChipEffects(this);
  c(getTags(ces), this$.tags);
})



###########################################################################/**
# @RdocMethod getFullName
#
# @title "Gets the full name of the output data set"
#
# \description{
#  @get "title", which is the name with comma separated tags.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFullName", "Crlmm", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})



###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path of this model"
#
# \description{
#  @get "title" where the parameter files are stored.
#  If non-existing, then the directory is created.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPath", "Crlmm", function(this, ...) {
  # Create the (sub-)directory tree for the dataset

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  ces <- getChipEffects(this);
  cdf <- getCdf(ces);
  chipType <- getChipType(cdf);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
})


###########################################################################/**
# @RdocMethod getChipEffects
#
# @title "Gets the chip effects for this model"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @see "SnpChipEffectSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChipEffects", "Crlmm", function(this, ...) {
  this$ces;
})


###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF structure for this model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCdfFile" object.
# }
#
# @author
#
# \seealso{
#   @seemethod "setCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCdf", "Crlmm", function(this, ...) {
  getCdf(this$ces);
})


###########################################################################/**
# @RdocMethod fit
#
# @title "Fits the model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{recalibrate}{If @TRUE, a second round of the fitting is done.}
#   \item{...}{Additional arguments.}
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
#
# @keyword IO
#*/###########################################################################
setMethodS3("fit", "Crlmm", function(this, recalibrate=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Requirements
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  require(oligo) || throw("Package 'oligo' not loaded.");

  # Assert that the PD enviroment package for the CDF is installed
  cdf <- getCdf(this);
  chipType <- getChipType(cdf);
  chipType <- gsub("-monocell$", "", chipType);
  pkgName <- sprintf("pd%s", cleancdfname(chipType, add=FALSE));
  if (!pkgName %in% .packages(all.available=TRUE)) {
    throw("Platform Design Environment (PDEnv) package '", 
                                              pkgName, "' not available.");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Default settings
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
  ces <- getChipEffects(this);
  qs <- extractSnpQSet(ces, transform=this$transform, verbose=less(verbose));
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
  key <- list(method="fitAffySnpMixture", path=getPath(ces), 
                                         sampleNames=getSampleNames(ces));
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
      bad <- which(myCalls == kk & llr < this$minLLRforCalls[kk]);
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
})


############################################################################
# HISTORY:
# 2006-11-20
# o Converted fitCrlmm() for SnpChipEffectSet into a model class.
# 2006-10-03
# o Created based on oligo::crlmm() v0.99.20.
############################################################################
