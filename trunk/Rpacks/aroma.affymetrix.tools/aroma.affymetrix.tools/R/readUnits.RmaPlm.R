setMethodS3("readUnits","RmaPlm",function(this, units=NULL, includeChipEffect=TRUE,includeProbeAffinity=TRUE,..., verbose=FALSE)
{
  verbose <- Arguments$getVerbose(verbose);

  ds <- getDataSet(this);
  verbose && enter(verbose, "Reading probe intensities from ", length(ds), " arrays");

  # Get the CDF cell indices
  verbose && enter(verbose, "Identifying CDF cell indices");
  cdfUnits <- getCellIndices(this, units=units, ...);
  verbose && print(verbose, cdfUnits[1]);
  verbose && exit(verbose);

  # Get the CEL intensities by units
  res <- getUnitIntensities(ds, units=cdfUnits, ...);
  verbose && str(verbose, res[1]);
  verbose && exit(verbose);
  
  if(!includeChipEffect){
    ces<- getChipEffectSet(this)
    cesData<-readUnits(ces,units=units)
    res<-mapply(res,cesData,SIMPLIFY=FALSE,FUN=function(x1,y1){
        out<-mapply(x1,y1,SIMPLIFY=FALSE,FUN=function(x2,y2){
            mat<-log2(x2$intensities); 
            chip<-log2(y2$theta); 
            mat<-sweep(mat,2,chip);
            mat<-2^mat
            return(list(intensities=mat))
        })
        names(out)<-names(x1)
        return(out)
    })
  }
  if(!includeProbeAffinity){
    paf<- getProbeAffinityFile(this)
    pafData<-readUnits(paf,units=units)
    res<-mapply(res,pafData,SIMPLIFY=FALSE,FUN=function(x1,y1){
        out<-mapply(x1,y1,SIMPLIFY=FALSE,FUN=function(x2,y2){
            mat<-log2(x2$intensities); 
            probe<-log2(y2$phi); 
            mat<-sweep(mat,1,probe);
            mat<-2^mat
            return(list(intensities=mat))})
        names(out)<-names(x1)
        return(out)
    })
  }
  res;

})
