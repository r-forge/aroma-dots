setMethodS3("mergeGroups","AffymetrixCelSet",function(this,units,arrays=NULL,fields=NULL,force=F,logscale=F,...){
    unitFcn<-getUnitMergeGroupsFunction(this,logscale=logscale,...) #a function that takes a groupList, arrays, and fields
    data<-readUnits(this,units=units,force=force)
    if(is.null(arrays)) arrays<-1:nbrOfArrays(this)
    if(!is.null(fields)) mergeData<-lapply(data,unitFcn,fields=fields,arrays=arrays)
    else mergeData<-lapply(data,unitFcn,arrays=arrays)
    return(mergeData)
})
